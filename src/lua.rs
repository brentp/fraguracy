use std::{collections::HashMap, ops::AddAssign};

use anyhow::Result;
use mlua::{prelude::LuaError, Function, Lua, UserData, UserDataFields, UserDataMethods, Value};
use rust_htslib::bam::{
    record::{Aux, Cigar},
    Record,
};

#[derive(Clone)]
pub struct LuaReadFilter {
    lua: Lua,
    filter_func: Function,
}

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

impl LuaReadFilter {
    pub fn skip_read(&self, read: &Record) -> Result<bool> {
        let r = self.lua.scope(|scope| {
            let globals = self.lua.globals();
            let user_data = scope.create_any_userdata_ref(read)?;
            globals.set("read", user_data).expect("failed to set read");
            self.filter_func.call::<bool>(())
        })?;
        Ok(r)
    }

    pub fn new(expression: &str, lua: Lua) -> Result<Self> {
        if !expression.contains("return") {
            return Err(anyhow::anyhow!(
                "expression must contain a return statement"
            ));
        }
        let filter_func = lua.load(expression).into_function()?;

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
                    Aux::Char(v) => Value::String(lua.create_string(&[v])?),
                    Aux::I8(v) => Value::Number(v as f64),
                    Aux::U8(v) => Value::Number(v as f64),
                    Aux::I16(v) => Value::Number(v as f64),
                    Aux::U16(v) => Value::Number(v as f64),
                    Aux::I32(v) => Value::Number(v as f64),
                    Aux::U32(v) => Value::Number(v as f64),
                    Aux::Float(v) => Value::Number(v as f64),
                    Aux::Double(v) => Value::Number(v as f64),
                    Aux::String(v) => Value::String(lua.create_string(&v)?),
                    Aux::ArrayFloat(v) => {
                        let mut arr = Vec::new();
                        for i in 0..v.len() {
                            arr.push(v.get(i).unwrap_or(f32::NAN) as f32);
                        }
                        Value::Table(lua.create_sequence_from(arr)?)
                    }
                    Aux::ArrayI32(v) => {
                        let mut arr = Vec::new();
                        for i in 0..v.len() {
                            arr.push(v.get(i).unwrap_or(i32::MIN) as i32);
                        }
                        Value::Table(lua.create_sequence_from(arr)?)
                    }
                    Aux::ArrayI8(v) => {
                        let mut arr = Vec::new();
                        for i in 0..v.len() {
                            arr.push(v.get(i).unwrap_or(i8::MIN) as i8);
                        }
                        Value::Table(lua.create_sequence_from(arr)?)
                    }
                    Aux::ArrayU8(v) => {
                        let mut arr = Vec::new();
                        for i in 0..v.len() {
                            arr.push(v.get(i).unwrap_or(u8::MIN) as u8);
                        }
                        Value::Table(lua.create_sequence_from(arr)?)
                    }
                    Aux::ArrayU16(v) => {
                        let mut arr = Vec::new();
                        for i in 0..v.len() {
                            arr.push(v.get(i).unwrap_or(u16::MIN) as u16);
                        }
                        Value::Table(lua.create_sequence_from(arr)?)
                    }
                    Aux::ArrayU32(v) => {
                        let mut arr = Vec::new();
                        for i in 0..v.len() {
                            arr.push(v.get(i).unwrap_or(u32::MIN) as u32);
                        }
                        Value::Table(lua.create_sequence_from(arr)?)
                    }
                    Aux::ArrayI16(v) => {
                        let mut arr = Vec::new();
                        for i in 0..v.len() {
                            arr.push(v.get(i).unwrap_or(i16::MIN) as i16);
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
            /*
            reg.add_field_function_get("bq", |_, this| {
                let qpos: usize = match this.named_user_value("qpos") {
                    Ok(qpos) => qpos,
                    Err(_) => {
                        return Ok(-1);
                    }
                };
                let this = this.borrow::<Record>()?;
                Ok(this.qual()[qpos] as i32)
            });
            reg.add_field_function_get("distance_from_5prime", |_, this| {
                let qpos: usize = match this.named_user_value("qpos") {
                    Ok(qpos) => qpos,
                    Err(_) => {
                        return Ok(-1);
                    }
                };
                let this = this.borrow::<Record>()?;
                if this.is_reverse() {
                    Ok(this.seq_len() as i32 - qpos as i32)
                } else {
                    Ok(qpos as i32)
                }
            });
            reg.add_field_function_get("distance_from_3prime", |_, this| {
                let qpos: usize = match this.named_user_value("qpos") {
                    Ok(qpos) => qpos,
                    Err(_) => {
                        return Ok(usize::MAX);
                    }
                };
                let this = this.borrow::<Record>()?;
                if this.is_reverse() {
                    Ok(qpos)
                } else {
                    Ok(this.seq_len() - qpos)
                }
            });
            */
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

        Ok(Self { lua, filter_func })
    }
}
