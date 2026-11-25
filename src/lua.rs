use std::{collections::HashMap, ops::AddAssign};

use anyhow::Result;
use mlua::{
    prelude::LuaError, AnyUserData, Compiler, Function, Lua, UserData, UserDataFields,
    UserDataMethods, Value,
};
use rust_htslib::bam::{
    record::{Aux, Cigar},
    Record,
};

pub struct LuaFilter {
    lua: Lua,
    read_func: Option<Function>,
    base_func: Option<Function>,
    chunk_func: Option<Function>,
    base_uses_read: bool,
    base_uses_base_char: bool,
}

#[derive(Debug, PartialEq, Eq)]
pub enum ChunkSkip {
    /// No bases in the chunk were skipped by the base Lua filter.
    None,
    /// Every eligible base in the chunk was skipped.
    All,
    /// Per-offset mask populated into the caller's scratch buffer (`1` = skip, `0` = keep).
    Mask,
}

pub struct ChunkEval<'a> {
    pub a: &'a Record,
    pub b: &'a Record,
    pub a_q0: u32,
    pub b_q0: u32,
    pub g0: u32,
    pub len: usize,
    pub min_bq: u8,
    pub a_quals: &'a [u8],
    pub b_quals: &'a [u8],
    pub a_seq: Option<&'a [u8]>,
    pub b_seq: Option<&'a [u8]>,
}

#[derive(Clone, Copy, Debug)]
pub struct BasePos {
    pub qpos: usize,
    pub gpos: u32,
    pub seq_len: usize,
    pub bq: u8,
    pub base: u8,
    pub is_reverse: bool,
    pub is_indel: bool,
    pub indel_length: i32,
}

pub const BASE_IN_READ_EXPR_ERR: &str =
    "read-level -l/--lua-expression cannot access `base.*`; use -L/--lua-base-expression instead";

pub struct Flags {
    flag: u16,
}

impl UserData for Flags {
    fn add_fields<M: UserDataFields<Self>>(fields: &mut M) {
        fields.add_field_method_get("paired", |_, this| Ok((this.flag & 0x1) != 0));
        fields.add_field_method_get("proper_pair", |_, this| Ok((this.flag & 0x2) != 0));
        fields.add_field_method_get("unmapped", |_, this| Ok((this.flag & 0x4) != 0));
        fields.add_field_method_get("mate_unmapped", |_, this| Ok((this.flag & 0x8) != 0));
        fields.add_field_method_get("reverse", |_, this| Ok((this.flag & 0x10) != 0));
        fields.add_field_method_get("forward", |_, this| Ok((this.flag & 0x10) == 0));
        fields.add_field_method_get("mate_reverse", |_, this| Ok((this.flag & 0x20) != 0));
        fields.add_field_method_get("mate_forward", |_, this| Ok((this.flag & 0x20) == 0));
        fields.add_field_method_get("read_1", |_, this| Ok((this.flag & 0x40) != 0));
        fields.add_field_method_get("read_2", |_, this| Ok((this.flag & 0x80) != 0));
        fields.add_field_method_get("secondary", |_, this| Ok((this.flag & 0x100) != 0));
        fields.add_field_method_get("primary", |_, this| Ok((this.flag & 0x100) == 0));
        fields.add_field_method_get("qcfail", |_, this| Ok((this.flag & 0x200) != 0));
        fields.add_field_method_get("duplicate", |_, this| Ok((this.flag & 0x400) != 0));
        fields.add_field_method_get("supplementary", |_, this| Ok((this.flag & 0x800) != 0));

        fields.add_field_method_get("flag", |_, this| Ok(this.flag));
    }
}

#[inline]
fn is_ident_char(b: u8) -> bool {
    b.is_ascii_alphanumeric() || b == b'_'
}

/// Returns true if `needle` appears as a standalone token.
/// Rejecting a preceding `.` prevents matching field names like `read.base_counts`
/// or `foo.read` as references to the top-level `base` or `read` arguments.
fn contains_token(expr: &str, needle: &str) -> bool {
    let bytes = expr.as_bytes();
    let nlen = needle.len();
    let mut start = 0;
    while let Some(pos) = expr[start..].find(needle) {
        let abs = start + pos;
        let before_ok = abs == 0 || (!is_ident_char(bytes[abs - 1]) && bytes[abs - 1] != b'.');
        let after = abs + nlen;
        let after_ok = after == bytes.len() || !is_ident_char(bytes[after]);
        if before_ok && after_ok {
            return true;
        }
        start = abs + 1;
    }
    false
}

/// Returns true if the `-L` expression may read `base.base` (either directly, via dynamic
/// indexing `base[...]`, via whitespace like `base .base`, or indirectly via aliasing `local b = base`
/// or passing `base` to a function). Only when every standalone `base` token is immediately followed
/// by a known non-`base` field can we safely skip decoding the sequence bytes.
fn base_needs_seq(expr: &str) -> bool {
    const NON_SEQ_FIELDS: [&str; 7] = [
        ".qpos",
        ".gpos",
        ".bq",
        ".distance_from_5prime",
        ".distance_from_3prime",
        ".is_indel",
        ".indel_length",
    ];
    let bytes = expr.as_bytes();
    let mut start = 0;
    while let Some(pos) = expr[start..].find("base") {
        let abs = start + pos;
        let before_ok = abs == 0 || (!is_ident_char(bytes[abs - 1]) && bytes[abs - 1] != b'.');
        let after = abs + 4;
        let after_ok = after == bytes.len() || !is_ident_char(bytes[after]);
        if before_ok && after_ok {
            let tail = &expr[after..];
            let tail_bytes = tail.as_bytes();
            let known_non_seq = NON_SEQ_FIELDS.iter().any(|f| {
                tail.starts_with(f)
                    && (tail_bytes.len() == f.len() || !is_ident_char(tail_bytes[f.len()]))
            });
            if !known_non_seq {
                return true;
            }
        }
        start = abs + 1;
    }
    false
}

/// Generate a Luau module returning `(eval_chunk, eval_single)` backed by the same `_user_pred`
/// and identical `base` table layout (so chunked bases and single-position indels cannot drift).
/// Always populates all numeric fields of `ba`/`bb` (so aliasing `local b = base`, passing `base`
/// to helpers, or `base .bq` always sees real values), and populates `.base` whenever
/// `uses_base_char` is true. Uses unaliased `buffer.readu8` / `string.sub` (required for Luau
/// `FASTCALL` opcodes under `lua.sandbox(true)` + opt-level 2) and a lazily cleared `mask_buf`.
///
/// Chunk return protocol: `nil` = 0 bases skipped, `true` = all eligible bases skipped,
/// `string` = `n`-byte mask (`1` = skip, `0` = keep).
fn build_base_lua(expr: &str, uses_base_char: bool) -> String {
    let (set_a_base, set_b_base) = if uses_base_char {
        (
            "            ba.base = string.sub(a_seq, off + 1, off + 1)\n",
            "                bb.base = string.sub(b_seq, off + 1, off + 1)\n",
        )
    } else {
        ("", "")
    };

    format!(
        r#"
local function _user_pred(read, base)
{expr}
end
local ba = {{ qpos = 0, gpos = 0, bq = 0, distance_from_5prime = 0, distance_from_3prime = 0, base = "", is_indel = false, indel_length = 0 }}
local bb = {{ qpos = 0, gpos = 0, bq = 0, distance_from_5prime = 0, distance_from_3prime = 0, base = "", is_indel = false, indel_length = 0 }}
local bs = {{ qpos = 0, gpos = 0, bq = 0, distance_from_5prime = 0, distance_from_3prime = 0, base = "", is_indel = false, indel_length = 0 }}
local mask_buf = buffer.create(1024)
local mask_cap = 1024
local function eval_chunk(read_a, read_b, a_q0, b_q0, g0, n, a_len, b_len, a_rev, b_rev, min_bq, a_qbuf, b_qbuf, a_seq, b_seq)
    if n > mask_cap then
        mask_buf = buffer.create(n)
        mask_cap = n
    end
    local a_max = a_len - 1
    local b_max = b_len - 1
    local a_d5_0 = if a_rev then a_max - a_q0 else a_q0
    local a_d5_step = if a_rev then -1 else 1
    local a_d3_0 = if a_rev then a_q0 else a_max - a_q0
    local a_d3_step = if a_rev then 1 else -1
    local b_d5_0 = if b_rev then b_max - b_q0 else b_q0
    local b_d5_step = if b_rev then -1 else 1
    local b_d3_0 = if b_rev then b_q0 else b_max - b_q0
    local b_d3_step = if b_rev then 1 else -1
    local skipped = 0
    local evaluated = 0
    local buf_cleared = false
    for off = 0, n - 1 do
        local aq = buffer.readu8(a_qbuf, off)
        local bq = buffer.readu8(b_qbuf, off)
        if aq >= min_bq and bq >= min_bq then
            evaluated += 1
            ba.qpos = a_q0 + off
            ba.gpos = g0 + off
            ba.bq = aq
            ba.distance_from_5prime = a_d5_0 + off * a_d5_step
            ba.distance_from_3prime = a_d3_0 + off * a_d3_step
{set_a_base}            local is_skip = _user_pred(read_a, ba)
            if not is_skip then
                bb.qpos = b_q0 + off
                bb.gpos = g0 + off
                bb.bq = bq
                bb.distance_from_5prime = b_d5_0 + off * b_d5_step
                bb.distance_from_3prime = b_d3_0 + off * b_d3_step
{set_b_base}                is_skip = _user_pred(read_b, bb)
            end
            if is_skip then
                if not buf_cleared then
                    buffer.fill(mask_buf, 0, 0, n)
                    buf_cleared = true
                end
                buffer.writeu8(mask_buf, off, 1)
                skipped += 1
            end
        end
    end
    if skipped == 0 then
        return nil
    elseif skipped == evaluated then
        return true
    else
        return buffer.readstring(mask_buf, 0, n)
    end
end
local function eval_single(read, qpos, gpos, seq_len, bq, base_str, is_rev, is_indel, indel_length)
    local max_q = if seq_len > 0 then seq_len - 1 else 0
    local tail = if max_q >= qpos then max_q - qpos else 0
    bs.qpos = qpos
    bs.gpos = gpos
    bs.bq = bq
    bs.distance_from_5prime = if is_rev then tail else qpos
    bs.distance_from_3prime = if is_rev then qpos else tail
    bs.base = base_str
    bs.is_indel = is_indel
    bs.indel_length = indel_length
    return _user_pred(read, bs)
end
return eval_chunk, eval_single
"#
    )
}

impl LuaFilter {
    #[inline(always)]
    pub fn has_read_filter(&self) -> bool {
        self.read_func.is_some()
    }

    #[inline(always)]
    pub fn has_base_filter(&self) -> bool {
        self.base_func.is_some()
    }

    #[cfg(test)]
    pub fn skip_read(&self, read: &Record) -> Result<bool> {
        let Some(ref func) = self.read_func else {
            return Ok(false);
        };
        let r = self.lua.scope(|scope| {
            let user_data = scope.create_any_userdata_ref(read)?;
            func.call::<bool>(user_data)
        })?;
        Ok(r)
    }

    pub fn skip_pair(&self, a: &Record, b: &Record) -> Result<bool> {
        let Some(ref func) = self.read_func else {
            return Ok(false);
        };
        let r = self.lua.scope(|scope| {
            let ud_a = scope.create_any_userdata_ref(a)?;
            if func.call::<bool>(ud_a)? {
                return Ok(true);
            }
            let ud_b = scope.create_any_userdata_ref(b)?;
            func.call::<bool>(ud_b)
        })?;
        Ok(r)
    }

    pub fn skip_base(&self, read: &Record, base: BasePos) -> Result<bool> {
        let Some(ref func) = self.base_func else {
            return Ok(false);
        };
        let base_bytes = if base.base != 0 {
            std::slice::from_ref(&base.base)
        } else {
            &[]
        };
        let base_str = self.lua.create_string(base_bytes)?;
        let invoke = |r: Value| -> mlua::Result<bool> {
            func.call::<bool>((
                r,
                base.qpos,
                base.gpos,
                base.seq_len,
                base.bq,
                base_str,
                base.is_reverse,
                base.is_indel,
                base.indel_length,
            ))
        };
        let r = if self.base_uses_read {
            self.lua.scope(|scope| {
                let ud_read = scope.create_any_userdata_ref(read)?;
                invoke(Value::UserData(ud_read))
            })?
        } else {
            invoke(Value::Nil)?
        };
        Ok(r)
    }

    /// Evaluate the per-base Lua expression across an entire overlapping chunk in one call.
    pub fn eval_chunk(&self, args: ChunkEval<'_>, out_mask: &mut Vec<u8>) -> Result<ChunkSkip> {
        let Some(ref func) = self.chunk_func else {
            return Ok(ChunkSkip::None);
        };
        if args.len == 0 {
            return Ok(ChunkSkip::None);
        }

        let a_quals_lua = self.lua.create_buffer(args.a_quals)?;
        let b_quals_lua = self.lua.create_buffer(args.b_quals)?;
        let a_seq_lua = match args.a_seq {
            Some(s) => Value::String(self.lua.create_string(s)?),
            None => Value::Nil,
        };
        let b_seq_lua = match args.b_seq {
            Some(s) => Value::String(self.lua.create_string(s)?),
            None => Value::Nil,
        };
        let a_len = args.a.seq_len();
        let b_len = args.b.seq_len();
        let a_rev = args.a.is_reverse();
        let b_rev = args.b.is_reverse();

        let invoke = |ra: Value, rb: Value| -> mlua::Result<Value> {
            func.call::<Value>((
                ra,
                rb,
                args.a_q0,
                args.b_q0,
                args.g0,
                args.len,
                a_len,
                b_len,
                a_rev,
                b_rev,
                args.min_bq,
                a_quals_lua,
                b_quals_lua,
                a_seq_lua,
                b_seq_lua,
            ))
        };

        let val: Value = if self.base_uses_read {
            self.lua.scope(|scope| {
                let ud_a: AnyUserData = scope.create_any_userdata_ref(args.a)?;
                let ud_b: AnyUserData = scope.create_any_userdata_ref(args.b)?;
                invoke(Value::UserData(ud_a), Value::UserData(ud_b))
            })?
        } else {
            invoke(Value::Nil, Value::Nil)?
        };

        match val {
            Value::Nil => Ok(ChunkSkip::None),
            Value::Boolean(true) => Ok(ChunkSkip::All),
            Value::String(s) => {
                out_mask.clear();
                out_mask.extend_from_slice(&s.as_bytes());
                Ok(ChunkSkip::Mask)
            }
            other => Err(anyhow::anyhow!(
                "unexpected return value from chunk base filter: {other:?}"
            )),
        }
    }

    #[inline(always)]
    pub fn needs_base_char(&self) -> bool {
        self.base_uses_base_char
    }

    pub fn new(
        read_expression: Option<&str>,
        base_expression: Option<&str>,
        lua: Lua,
    ) -> Result<Self> {
        // Enable Luau safeenv sandbox + opt-level 2 so globals/builtins (`buffer`, `string`,
        // `math`) are read-only and compile to native FASTCALL VM opcodes instead of table lookups.
        lua.sandbox(true)?;
        lua.set_compiler(Compiler::new().set_optimization_level(2));

        let read_func = if let Some(expr) = read_expression {
            if !expr.contains("return") {
                return Err(anyhow::anyhow!(
                    "lua expression must contain a return statement"
                ));
            }
            // Sentinel `base` upvalue gives a clear runtime error if a `-l` expression accesses
            // `base.*`, while any user `local base = ...` naturally shadows it via Lua scoping.
            let wrapped = format!(
                "local base = setmetatable({{}}, {{ __index = function() error(\"{BASE_IN_READ_EXPR_ERR}\") end }})\nreturn function(read)\n{expr}\nend"
            );
            Some(lua.load(&wrapped).eval::<Function>()?)
        } else {
            None
        };

        let mut base_uses_read = false;
        let mut base_uses_base_char = false;
        let (base_func, chunk_func) = if let Some(expr) = base_expression {
            if !expr.contains("return") {
                return Err(anyhow::anyhow!(
                    "lua base expression must contain a return statement"
                ));
            }
            base_uses_read = contains_token(expr, "read");
            base_uses_base_char = base_needs_seq(expr);

            let (chunk, single): (Function, Function) = lua
                .load(build_base_lua(expr, base_uses_base_char))
                .eval()?;
            (Some(single), Some(chunk))
        } else {
            (None, None)
        };

        lua.register_userdata_type::<Record>(|reg| {
            reg.add_field_method_get("mapping_quality", |_, this| Ok(this.mapq()));
            reg.add_field_method_get("flags", |_, this| Ok(Flags { flag: this.flags() }));
            reg.add_field_method_get("tid", |_, this| Ok(this.tid()));
            reg.add_field_method_get("start", |_, this| Ok(this.pos()));
            reg.add_field_method_get("stop", |_, this| Ok(this.cigar().end_pos()));
            reg.add_field_method_get("length", |_, this| Ok(this.seq_len()));
            reg.add_field_method_get("insert_size", |_, this| Ok(this.insert_size()));
            reg.add_field_method_get("qname", |_, this| {
                let q = this.qname();
                Ok(std::str::from_utf8(q).unwrap_or("").to_string())
            });
            reg.add_field_method_get("sequence", |_, this| {
                let seq = this.seq();
                Ok(std::str::from_utf8(&seq.as_bytes())
                    .unwrap_or("")
                    .to_string())
            });

            reg.add_field_method_get("soft_clips_3_prime", |_, this| {
                let cigar = this.cigar();
                if this.is_reverse() {
                    Ok(cigar.leading_softclips())
                } else {
                    Ok(cigar.trailing_softclips())
                }
            });
            reg.add_field_method_get("soft_clips_5_prime", |_, this| {
                let cigar = this.cigar();
                if this.is_reverse() {
                    Ok(cigar.trailing_softclips())
                } else {
                    Ok(cigar.leading_softclips())
                }
            });
            reg.add_field_method_get("average_base_quality", |_, this| {
                let qual = this.qual();
                let sum = qual.iter().map(|q| *q as u64).sum::<u64>();
                let count = qual.len();
                Ok(sum as f64 / count as f64)
            });

            reg.add_method("tag", |lua, this: &Record, tag: String| {
                let tag = tag.as_bytes();
                let aux = this.aux(tag).map_err(LuaError::external)?;
                let lua_val: Value = match aux {
                    Aux::Char(v) => Value::String(lua.create_string([v])?),
                    Aux::I8(v) => Value::Number(v as f64),
                    Aux::U8(v) => Value::Number(v as f64),
                    Aux::I16(v) => Value::Number(v as f64),
                    Aux::U16(v) => Value::Number(v as f64),
                    Aux::I32(v) => Value::Number(v as f64),
                    Aux::U32(v) => Value::Number(v as f64),
                    Aux::Float(v) => Value::Number(v as f64),
                    Aux::Double(v) => Value::Number(v),
                    Aux::String(v) => Value::String(lua.create_string(v)?),
                    Aux::ArrayFloat(v) => {
                        let mut arr = Vec::new();
                        for i in 0..v.len() {
                            arr.push(v.get(i).unwrap_or(f32::NAN));
                        }
                        Value::Table(lua.create_sequence_from(arr)?)
                    }
                    Aux::ArrayI32(v) => {
                        let mut arr = Vec::new();
                        for i in 0..v.len() {
                            arr.push(v.get(i).unwrap_or(i32::MIN));
                        }
                        Value::Table(lua.create_sequence_from(arr)?)
                    }
                    Aux::ArrayI8(v) => {
                        let mut arr = Vec::new();
                        for i in 0..v.len() {
                            arr.push(v.get(i).unwrap_or(i8::MIN));
                        }
                        Value::Table(lua.create_sequence_from(arr)?)
                    }
                    Aux::ArrayU8(v) => {
                        let mut arr = Vec::new();
                        for i in 0..v.len() {
                            arr.push(v.get(i).unwrap_or(u8::MIN));
                        }
                        Value::Table(lua.create_sequence_from(arr)?)
                    }
                    Aux::ArrayU16(v) => {
                        let mut arr = Vec::new();
                        for i in 0..v.len() {
                            arr.push(v.get(i).unwrap_or(u16::MIN));
                        }
                        Value::Table(lua.create_sequence_from(arr)?)
                    }
                    Aux::ArrayU32(v) => {
                        let mut arr = Vec::new();
                        for i in 0..v.len() {
                            arr.push(v.get(i).unwrap_or(u32::MIN));
                        }
                        Value::Table(lua.create_sequence_from(arr)?)
                    }
                    Aux::ArrayI16(v) => {
                        let mut arr = Vec::new();
                        for i in 0..v.len() {
                            arr.push(v.get(i).unwrap_or(i16::MIN));
                        }
                        Value::Table(lua.create_sequence_from(arr)?)
                    }
                    Aux::HexByteArray(v) => {
                        let lstr = String::from_utf8_lossy(v.as_bytes()).to_string();
                        Value::String(lua.create_string(&lstr)?)
                    }
                };
                Ok(Some(lua_val))
            });

            // count the number of A, C, G, T, N in the read. Always capitalize and return a table
            reg.add_field_method_get("base_counts", |_, this| {
                let seq = this.seq();
                let mut counts = HashMap::new();
                for i in 0..seq.len() {
                    let base = seq[i].to_ascii_uppercase();
                    counts.entry(base).or_insert(0).add_assign(1);
                }
                Ok(counts)
            });
            reg.add_field_method_get("n_proportion", |_, this| {
                let seq = this.seq();
                let mut count = 0;
                for i in 0..seq.len() {
                    let base = seq[i].to_ascii_uppercase();
                    if base == b'N' {
                        count += 1;
                    }
                }
                Ok(count as f64 / seq.len() as f64)
            });

            reg.add_method("n_proportion_3_prime", |_, this, n_bases: usize| {
                let seq = this.seq();
                let mut count = 0;
                let reverse = this.is_reverse();
                for i in 0..n_bases {
                    let base =
                        seq[if reverse { i } else { seq.len() - 1 - i }].to_ascii_uppercase();
                    if base == b'N' {
                        count += 1;
                    }
                }
                Ok(count as f64 / n_bases as f64)
            });

            reg.add_method("n_proportion_5_prime", |_, this, n_bases: usize| {
                let seq = this.seq();
                let mut count = 0;
                let reverse = this.is_reverse();
                for i in 0..n_bases {
                    let base =
                        seq[if reverse { seq.len() - 1 - i } else { i }].to_ascii_uppercase();
                    if base == b'N' {
                        count += 1;
                    }
                }
                Ok(count as f64 / n_bases as f64)
            });

            reg.add_field_method_get("indel_count", |_, this| {
                let cigar = this.cigar();
                let mut count = 0;
                for op in cigar.iter() {
                    count += match op {
                        Cigar::Ins(_len) => 1,
                        Cigar::Del(_len) => 1,
                        _ => 0,
                    }
                }
                Ok(count)
            });
        })?;

        Ok(Self {
            lua,
            read_func,
            base_func,
            chunk_func,
            base_uses_read,
            base_uses_base_char,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{Cigar, CigarString};

    fn make_record(qname: &[u8], seq: &[u8], qual: &[u8], mapq: u8, reverse: bool) -> Record {
        let mut rec = Record::new();
        rec.set(
            qname,
            Some(&CigarString(vec![Cigar::Match(seq.len() as u32)])),
            seq,
            qual,
        );
        rec.set_mapq(mapq);
        if reverse {
            rec.set_reverse();
        }
        rec
    }

    #[test]
    fn test_lua_filter_local_arg_and_pair() {
        let filter = LuaFilter::new(
            Some(
                "-- check mapq\nlocal mq = read.mapping_quality\nreturn mq < 30 -- trailing comment",
            ),
            None,
            Lua::new(),
        )
        .unwrap();

        let high = make_record(b"r1", b"ACGTACGT", &[30; 8], 60, false);
        let low = make_record(b"r2", b"ACGTACGT", &[30; 8], 10, false);

        assert!(!filter.skip_read(&high).unwrap());
        assert!(filter.skip_read(&low).unwrap());
        assert!(!filter.skip_pair(&high, &high).unwrap());
        assert!(filter.skip_pair(&low, &high).unwrap());
        assert!(filter.skip_pair(&high, &low).unwrap());
    }

    #[test]
    fn test_lua_filter_requires_return() {
        assert!(LuaFilter::new(Some("read.mapping_quality < 30"), None, Lua::new()).is_err());
        assert!(LuaFilter::new(None, Some("base.qpos < 5"), Lua::new()).is_err());
    }

    #[test]
    fn test_read_expression_rejects_base_access() {
        let r = make_record(b"r1", b"ACGTACGT", &[30; 8], 60, false);

        // Accessing `base.*` in `-l` fails with a clear sentinel error on evaluation.
        let bad = LuaFilter::new(Some("return base.qpos < 5"), None, Lua::new()).unwrap();
        let err = bad.skip_read(&r).unwrap_err().to_string();
        assert!(err.contains(BASE_IN_READ_EXPR_ERR), "got: {err}");

        let ok = LuaFilter::new(
            Some("return read.average_base_quality < 20 and read.base_counts ~= nil"),
            None,
            Lua::new(),
        )
        .unwrap();
        assert!(!ok.skip_read(&r).unwrap());

        // Local variable named `base` (including multi-binding and function form) shadows the sentinel.
        for expr in [
            "local base = read.base_counts; return base ~= nil",
            "local x, base = 1, read.base_counts; return (base.N or 0) > x",
            "local x: number, base = 1, read.base_counts; return (base.N or 0) > x",
            "local function base() return false end; return base()",
        ] {
            let f = LuaFilter::new(Some(expr), None, Lua::new()).unwrap();
            assert!(
                f.skip_read(&r).is_ok(),
                "expected `-l` to allow `{expr}`"
            );
        }
    }

    #[test]
    fn test_base_filter_single_and_chunk_equivalence() {
        for expr in [
            "return base.distance_from_3prime < 2 or base.bq < 20 or base.base == 'G'",
            "local d = base.distance_from_3prime\nreturn d < 2 or base.bq < 20 or base.base == 'G'",
            // Indirect `base` access: aliasing, helper function, and whitespace before `.`
            "local b = base; return b.distance_from_3prime < 2 or b.bq < 20 or b.base == 'G'",
            "local function f(x) return x.distance_from_3prime < 2 or x.bq < 20 or x.base == 'G' end\nreturn f(base)",
            "return base .distance_from_3prime < 2 or base .bq < 20 or base .base == 'G'",
            // Truthy non-boolean return value (must match Lua truthiness in both chunk and single/indel paths)
            "return (base.distance_from_3prime < 2 or base.bq < 20 or base.base == 'G') and 'skip'",
        ] {
            let filter = LuaFilter::new(None, Some(expr), Lua::new()).unwrap();
            assert!(!filter.base_uses_read);
            assert!(filter.needs_base_char());

            let a_seq = b"ACGTACGT";
            let b_seq = b"ACTTACCT";
            let a_qual = [30u8, 30, 30, 30, 15, 30, 30, 30];
            let b_qual = [30u8, 30, 30, 30, 30, 30, 30, 30];
            let a = make_record(b"r1", a_seq, &a_qual, 60, false);
            let b = make_record(b"r2", b_seq, &b_qual, 60, true);

            let mut mask = Vec::new();
            let res = filter
                .eval_chunk(
                    ChunkEval {
                        a: &a,
                        b: &b,
                        a_q0: 0,
                        b_q0: 0,
                        g0: 100,
                        len: 8,
                        min_bq: 5,
                        a_quals: &a_qual,
                        b_quals: &b_qual,
                        a_seq: Some(a_seq),
                        b_seq: Some(b_seq),
                    },
                    &mut mask,
                )
                .unwrap();
            assert_eq!(res, ChunkSkip::Mask);
            assert_eq!(mask, vec![1, 1, 1, 0, 1, 0, 1, 1]);

            for i in 0..8 {
                let skip_a = filter
                    .skip_base(
                        &a,
                        BasePos {
                            qpos: i,
                            gpos: 100 + i as u32,
                            seq_len: 8,
                            bq: a_qual[i],
                            base: a_seq[i],
                            is_reverse: false,
                            is_indel: false,
                            indel_length: 0,
                        },
                    )
                    .unwrap();
                let skip_b = filter
                    .skip_base(
                        &b,
                        BasePos {
                            qpos: i,
                            gpos: 100 + i as u32,
                            seq_len: 8,
                            bq: b_qual[i],
                            base: b_seq[i],
                            is_reverse: true,
                            is_indel: false,
                            indel_length: 0,
                        },
                    )
                    .unwrap();
                assert_eq!(mask[i] != 0, skip_a || skip_b);
            }
        }
    }

    #[test]
    fn test_base_filter_none_and_all_fast_paths() {
        let filter_none =
            LuaFilter::new(None, Some("return base.qpos > 1000"), Lua::new()).unwrap();
        let filter_all =
            LuaFilter::new(None, Some("return base.qpos >= 0"), Lua::new()).unwrap();
        let a = make_record(b"r1", b"ACGT", &[30; 4], 60, false);
        let b = make_record(b"r2", b"ACGT", &[30; 4], 60, true);
        let mut mask = Vec::new();
        let args = ChunkEval {
            a: &a,
            b: &b,
            a_q0: 0,
            b_q0: 0,
            g0: 10,
            len: 4,
            min_bq: 5,
            a_quals: &[30; 4],
            b_quals: &[30; 4],
            a_seq: None,
            b_seq: None,
        };

        assert_eq!(filter_none.eval_chunk(args, &mut mask).unwrap(), ChunkSkip::None);
        let args2 = ChunkEval {
            a: &a,
            b: &b,
            a_q0: 0,
            b_q0: 0,
            g0: 10,
            len: 4,
            min_bq: 5,
            a_quals: &[30; 4],
            b_quals: &[30; 4],
            a_seq: None,
            b_seq: None,
        };
        assert_eq!(filter_all.eval_chunk(args2, &mut mask).unwrap(), ChunkSkip::All);
    }

    #[test]
    fn test_base_filter_with_read_and_indel() {
        let filter = LuaFilter::new(
            Some("return read.mapping_quality < 30"),
            Some("return (read.flags.reverse and base.qpos < 2) or (base.is_indel and base.indel_length < -1)"),
            Lua::new(),
        )
        .unwrap();
        assert!(filter.has_read_filter());
        assert!(filter.has_base_filter());
        assert!(filter.base_uses_read);
        assert!(!filter.needs_base_char());

        let a = make_record(b"r1", b"ACGTAC", &[30; 6], 60, false);
        let b = make_record(b"r2", b"ACGTAC", &[30; 6], 60, true);
        assert!(!filter.skip_pair(&a, &b).unwrap());

        let mut mask = Vec::new();
        let res = filter
            .eval_chunk(
                ChunkEval {
                    a: &a,
                    b: &b,
                    a_q0: 0,
                    b_q0: 0,
                    g0: 50,
                    len: 6,
                    min_bq: 5,
                    a_quals: &[30; 6],
                    b_quals: &[30; 6],
                    a_seq: None,
                    b_seq: None,
                },
                &mut mask,
            )
            .unwrap();
        assert_eq!(res, ChunkSkip::Mask);
        assert_eq!(mask, vec![1, 1, 0, 0, 0, 0]);

        let indel_del2 = BasePos {
            qpos: 3,
            gpos: 53,
            seq_len: 6,
            bq: 30,
            base: b'A',
            is_reverse: false,
            is_indel: true,
            indel_length: -2,
        };
        assert!(filter.skip_base(&a, indel_del2).unwrap());
        assert!(!filter.skip_base(&a, BasePos { indel_length: 1, ..indel_del2 }).unwrap());
    }
}
