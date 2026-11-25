use bpci::*;
use ndarray::prelude::Array;
use ndarray::Array6;
use rust_htslib::bam::{
    record::{Cigar, CigarStringView},
    IndexedReader, Record,
};
use rust_lapper::Lapper;
use rustc_hash::FxHashMap;

use crate::depth::DepthTracker;
use crate::homopolymer as hp;
use crate::lua::{BasePos, ChunkEval, ChunkSkip, LuaFilter, BASE_IN_READ_EXPR_ERR};
use crate::pileup::PileupWindow;
use std::collections::HashMap;
use std::fmt;
use std::hash::Hash;

#[derive(Eq, Hash, PartialEq, Ord, PartialOrd)]
pub(crate) struct Position {
    pub tid: u16,
    pub pos: u32,
    pub bq_bin: u8,
}

type Length = i32;

pub(crate) const MAX_HP_DIST: i16 = 15;
const HP_BINS: usize = (2 * MAX_HP_DIST + 2) as usize; // 32
const NO_HP_BIN: usize = HP_BINS - 1; // 31 = "NA"
const CNTS_BIN_STRIDE: usize = 5 * 2 * HP_BINS; // 320

/// BAM 4-bit nucleotide encoding (`=ACMGRSVTWYHKDBN`) -> collapsed 2-base context:
/// A(1)/T(8) -> 0, C(2)/G(4) -> 1, anything else -> 2.
/// Lets the hot matching-base path avoid ASCII decoding + match branching.
const ENC_TO_CTX2: [usize; 16] = [2, 0, 1, 2, 1, 2, 2, 2, 0, 2, 2, 2, 2, 2, 2, 2];
const ENC_TO_BASE: [u8; 16] = *b"=ACMGRSVTWYHKDBN";

/// Flat C-order index into `InnerCounts::cnts` with shape
/// `(read12: 2, rev: 2, bin: bins, bq: 5, ctx2: 2, hp_dist: 32)`:
/// `row_base` encodes `(read12, rev)`, `bin_off` encodes `bin * 320`,
/// `bq << 6` encodes `bq * 64`, and `ctx2 << 5` encodes `ctx2 * 32`.
#[inline(always)]
fn cnts_flat_idx(row_base: usize, bin_off: usize, bq: u8, ctx2: usize, hp_dist: usize) -> usize {
    debug_assert!(bq < 5 && ctx2 < 2 && hp_dist < HP_BINS);
    row_base + bin_off + ((bq as usize) << 6) + (ctx2 << 5) + hp_dist
}

/// Map `(true_base, error_base)` to one of the 6 canonical pyrimidine/purine-collapsed
/// substitution contexts (or 6 for `NN`), avoiding `HashMap` lookups on mismatches.
#[inline(always)]
pub(crate) fn context_lookup(a: u8, b: u8) -> usize {
    match (a, b) {
        (b'T', b'G') | (b'A', b'C') => 0,
        (b'T', b'C') | (b'A', b'G') => 1,
        (b'T', b'A') | (b'A', b'T') => 2,
        (b'C', b'A') | (b'G', b'T') => 3,
        (b'C', b'G') | (b'G', b'C') => 4,
        (b'C', b'T') | (b'G', b'A') => 5,
        (b'N', b'N') => 6,
        _ => panic!("unexpected base context pair: ({}, {})", a as char, b as char),
    }
}

/// Fast check to skip `indel_error_pieces_with_qpos` allocations for ungapped reads.
#[inline(always)]
fn has_indel(cigar: &CigarStringView) -> bool {
    if cigar.len() == 1 && matches!(cigar[0], Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_)) {
        return false;
    }
    cigar
        .iter()
        .any(|op| matches!(op, Cigar::Ins(_) | Cigar::Del(_)))
}

/// Returns the homopolymer distance with the minimum absolute value.
/// If both distances are available, returns the one with smaller absolute value.
/// If only one is available, returns that one.
/// If neither is available, returns None.
fn min_abs_hp_distance(dist_a: Option<i16>, dist_b: Option<i16>) -> Option<i16> {
    match (dist_a, dist_b) {
        (Some(a), Some(b)) => {
            if a.abs() <= b.abs() {
                Some(a)
            } else {
                Some(b)
            }
        }
        (Some(a), None) => Some(a),
        (None, Some(b)) => Some(b),
        (None, None) => None,
    }
}

/// Result of checking a genomic overlap chunk `[g0, g1)` once against an include or exclude tree.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ChunkFilter {
    /// Every base in the chunk passes (no per-base check needed).
    PassAll,
    /// Every base in the chunk is filtered out (skip the entire chunk).
    SkipAll,
    /// Only part of the chunk overlaps; walk the merged sorted intervals with a cursor.
    Partial,
}

/// Chunk-level include/exclude region filter backed by a reusable sorted, merged interval buffer
/// and a monotonic per-base cursor.
struct RegionFilter {
    mode: ChunkFilter,
    ivs: Vec<(u32, u32)>,
    cursor: usize,
    include: bool,
}

impl RegionFilter {
    fn new(include: bool) -> Self {
        Self {
            mode: ChunkFilter::PassAll,
            ivs: Vec::new(),
            cursor: 0,
            include,
        }
    }

    /// Query `tree` once for chunk `[g0, g1)` and cache merged intervals if partially overlapping.
    fn prepare(&mut self, tree: Option<&Lapper<u32, u8>>, g0: u32, g1: u32) -> ChunkFilter {
        self.cursor = 0;
        let Some(t) = tree else {
            self.mode = ChunkFilter::PassAll;
            return self.mode;
        };
        let (on_empty, on_full) = if self.include {
            (ChunkFilter::SkipAll, ChunkFilter::PassAll)
        } else {
            (ChunkFilter::PassAll, ChunkFilter::SkipAll)
        };
        let mut iter = t.find(g0, g1);
        let Some(first) = iter.next() else {
            self.mode = on_empty;
            return self.mode;
        };
        if first.start <= g0 && first.stop >= g1 {
            self.mode = on_full;
            return self.mode;
        }
        self.ivs.clear();
        self.ivs.push((first.start, first.stop));
        for iv in iter {
            // Merge overlapping/bookended intervals so the per-base cursor advances monotonically.
            if let Some(last) = self.ivs.last_mut() {
                if iv.start <= last.1 {
                    last.1 = last.1.max(iv.stop);
                    continue;
                }
            }
            self.ivs.push((iv.start, iv.stop));
        }
        self.mode = if self.ivs.len() == 1 && self.ivs[0].0 <= g0 && self.ivs[0].1 >= g1 {
            on_full
        } else {
            ChunkFilter::Partial
        };
        self.mode
    }

    /// Return whether monotonically increasing genomic coordinate `pos` passes this filter.
    #[inline(always)]
    fn passes(&mut self, pos: u32) -> bool {
        if self.mode != ChunkFilter::Partial {
            return true;
        }
        while self.cursor < self.ivs.len() && self.ivs[self.cursor].1 <= pos {
            self.cursor += 1;
        }
        let in_iv = self.cursor < self.ivs.len() && pos >= self.ivs[self.cursor].0;
        in_iv == self.include
    }
}

/// Per-chromosome state updated once on chromosome boundaries in `process_bam`.
pub(crate) struct ChromCtx<'a> {
    pub(crate) chrom: &'a str,
    pub(crate) chrom_seq: Option<Vec<u8>>,
    pub(crate) include_tree: Option<&'a Lapper<u32, u8>>,
    pub(crate) exclude_tree: Option<&'a Lapper<u32, u8>>,
    pub(crate) hp_tree: Option<Lapper<u32, u8>>,
}

struct Mate {
    start: u32,
    end: u32,
    strand: i8,
    read12: usize,
    rev_idx: usize,
    row_base: usize,
}

impl Mate {
    #[inline(always)]
    fn new(r: &Record, cigar: &CigarStringView, bins: usize) -> Self {
        let is_rev = r.is_reverse();
        let read12 = 1 - r.is_first_in_template() as usize;
        let rev_idx = is_rev as usize;
        Self {
            start: r.pos() as u32,
            end: cigar.end_pos() as u32,
            strand: if is_rev { -1 } else { 1 },
            read12,
            rev_idx,
            row_base: (read12 * 2 + rev_idx) * bins * CNTS_BIN_STRIDE,
        }
    }
}

pub(crate) struct Counts {
    pub(crate) pileup: Option<PileupWindow>,
    //  read, f/r pos, bq, bp, ctx{6} */
    pub(crate) counts: InnerCounts,
    pub(crate) depth: Option<DepthTracker>,
    min_base_qual: u8,
    min_map_qual: u8,
    bin_size: u32,
    lua_filter: Option<LuaFilter>,
    logged_read_err: bool,
    logged_base_err: bool,
    /// Precomputed `(qpos / bin_size) * CNTS_BIN_STRIDE` lookup table (avoids hardware `div` per base).
    bin_off_table: Vec<usize>,
    base_skip_mask: Vec<u8>,
    a_seq_scratch: Vec<u8>,
    b_seq_scratch: Vec<u8>,
    pieces_scratch: Vec<[Coordinates; 3]>,
    inc_filter: RegionFilter,
    exc_filter: RegionFilter,
}

pub(crate) struct InnerCounts {
    // genome_pos
    pub(crate) errs: Array6<u64>,
    //  read, f/r, pos, bq, ctx{2}, hp_dist */
    pub(crate) cnts: Array6<u64>,
    pub(crate) mismatches: u64,
    pub(crate) matches: u64,

    // position -> error count. nice to find sites that are error-prone.
    pub(crate) error_positions: FxHashMap<Position, [u32; 7]>,
    // position -> indel error counts
    pub(crate) indel_error_positions: FxHashMap<(Position, Length, i16), u32>,
}

fn argmax<T: Ord>(slice: &[T]) -> Option<usize> {
    (0..slice.len()).max_by_key(|i| &slice[*i])
}

impl std::ops::AddAssign<InnerCounts> for InnerCounts {
    fn add_assign(&mut self, o: InnerCounts) {
        self.errs.add_assign(&o.errs);
        self.cnts.add_assign(&o.cnts);
        self.mismatches += o.mismatches;
        self.matches += o.matches;

        for (pos, cnt) in o.error_positions.into_iter() {
            let entry = self.error_positions.entry(pos).or_insert([0; 7]);
            for i in 0..entry.len() {
                entry[i] += cnt[i];
            }
        }
        for (pos, cnt) in o.indel_error_positions.into_iter() {
            *(self.indel_error_positions.entry(pos)).or_insert(0) += cnt;
        }
    }
}

pub(crate) struct Stat {
    pub ci: ConfidenceInterval,
    read12: u8,
    fr: u8,
    bq_bin: u8,
    read_pos: u32,
    context: [char; 2],
    homopolymer_distance: i16,
    total_count: u64,
    error_count: u64,
}

unsafe impl std::marker::Sync for Counts {}

impl fmt::Display for Stat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (lo, hi) = self.confidence_interval(&self.ci);
        let hp_dist_str = if self.homopolymer_distance == MAX_HP_DIST + 1 {
            "NA".to_string()
        } else {
            self.homopolymer_distance.to_string()
        };
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}{}\t{}\t{}\t{}\t{:e}\t{:e}",
            ["r1", "r2"][self.read12 as usize],
            ["f", "r"][self.fr as usize],
            Q_LOOKUP[self.bq_bin as usize],
            self.read_pos,
            self.context[0],
            self.context[1],
            hp_dist_str,
            self.total_count,
            self.error_count,
            lo.max(0.0),
            hi.max(0.0),
        )
    }
}

#[derive(Debug, Clone, clap::ValueEnum, Default)]
pub enum ConfidenceInterval {
    #[default]
    AgrestiCoull,
    Wald,
    Wilson,
    //WilsonWithCC,
}

impl fmt::Display for ConfidenceInterval {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl Stat {
    pub(crate) fn header() -> String {
        String::from(
            "read12\tFR\tbq_bin\tread_pos\tcontext\thp_dist\ttotal_count\terror_count\terr_rate_lo\terr_rate_hi",
        )
    }

    pub(crate) fn confidence_interval(&self, ci: &ConfidenceInterval) -> (f64, f64) {
        let sample = bpci::NSuccessesSample::new(self.total_count as f64, self.error_count as f64)
            .expect("error with proportion");

        let f = match ci {
            ConfidenceInterval::AgrestiCoull => sample.agresti_coull(1.960),
            ConfidenceInterval::Wald => sample.wald(1.960),
            ConfidenceInterval::Wilson => sample.wilson_score(1.960),
        };

        (f.lower(), f.upper())
    }

    pub(crate) fn from_counts(
        c: &InnerCounts,
        bin_size: usize,
        ci: ConfidenceInterval,
    ) -> Vec<Stat> {
        let mut stats = vec![];
        for readi in 0..c.cnts.shape()[0] {
            for fri in 0..c.cnts.shape()[1] {
                for read_posi in 0..c.cnts.shape()[2] {
                    for bqi in 0..c.cnts.shape()[3] {
                        for ctx6i in 0..c.errs.shape()[4] {
                            for hp_dist in 0..c.errs.shape()[5] {
                                let n_err = c.errs[[readi, fri, read_posi, bqi, ctx6i, hp_dist]];

                                // from ctx6i, we get the original context.
                                let bases = CONTEXT_TO_CONTEXT2[ctx6i];

                                let ctx2i = Counts::base_to_ctx2(bases[0] as u8);
                                let n_tot = c.cnts[[readi, fri, read_posi, bqi, ctx2i, hp_dist]];
                                if n_tot < n_err {
                                    eprintln!(
                                        "BAD: {ctx6i} -> {bases:?}. ctx2i:{ctx2i}",
                                        ctx6i = ctx6i,
                                        bases = bases,
                                        ctx2i = ctx2i
                                    );
                                }

                                stats.push(Stat {
                                    ci: ci.clone(),
                                    read12: readi as u8,
                                    fr: fri as u8,
                                    bq_bin: bqi as u8,
                                    read_pos: (read_posi * bin_size) as u32,
                                    context: bases,
                                    total_count: n_tot,
                                    error_count: n_err,
                                    homopolymer_distance: hp_dist as i16 - MAX_HP_DIST,
                                })
                            }
                        }
                    }
                }
            }
        }
        stats
    }
}

impl InnerCounts {
    pub(crate) fn new(bins: usize) -> Self {
        InnerCounts {
            cnts: Array::zeros((2, 2, bins, 5, 2, (2 * MAX_HP_DIST + 2) as usize)),
            errs: Array::zeros((2, 2, bins, 5, 6, (2 * MAX_HP_DIST + 2) as usize)),
            mismatches: 0,
            matches: 0,
            error_positions: FxHashMap::default(),
            indel_error_positions: FxHashMap::default(),
        }
    }
}

impl Counts {
    pub(crate) fn new(
        ir: Option<IndexedReader>,
        bins: usize,
        bin_size: u32,
        min_base_qual: u8,
        min_map_qual: u8,
        lua_filter: Option<LuaFilter>,
    ) -> Self {
        let bin_off_table = (0..bins * bin_size as usize)
            .map(|qpos| (qpos / bin_size as usize) * CNTS_BIN_STRIDE)
            .collect();
        Counts {
            pileup: ir.map(PileupWindow::new),
            counts: InnerCounts::new(bins),
            depth: None,
            min_base_qual,
            min_map_qual,
            bin_size,
            lua_filter,
            logged_read_err: false,
            logged_base_err: false,
            bin_off_table,
            base_skip_mask: Vec::new(),
            a_seq_scratch: Vec::new(),
            b_seq_scratch: Vec::new(),
            pieces_scratch: Vec::with_capacity(4),
            inc_filter: RegionFilter::new(true),
            exc_filter: RegionFilter::new(false),
        }
    }

    pub(crate) fn skip_pair(&mut self, a: &Record, b: &Record) -> bool {
        let Some(ref lf) = self.lua_filter else {
            return false;
        };
        if !lf.has_read_filter() {
            return false;
        }
        match lf.skip_pair(a, b) {
            Ok(s) => s,
            Err(e) => {
                let msg = e.to_string();
                if msg.contains(BASE_IN_READ_EXPR_ERR) {
                    panic!("{msg}");
                }
                if !self.logged_read_err {
                    log::error!(
                        "error evaluating user expression for read (further errors suppressed): {e}"
                    );
                    self.logged_read_err = true;
                }
                false
            }
        }
    }

    pub(crate) fn set_depth_writer(&mut self, path: &str) -> std::io::Result<()> {
        self.depth = Some(DepthTracker::new(path)?);
        Ok(())
    }

    #[inline(always)]
    fn qual_to_bin(q: u8) -> u8 {
        match q {
            0..=5 => 0,
            6..=19 => 1,
            20..=36 => 2,
            37..=59 => 3,
            _ => 4,
        }
    }

    #[inline(always)]
    fn base_to_ctx2(b: u8) -> usize {
        match b as char {
            'A' | 'T' => 0,
            'C' | 'G' => 1,
            'N' => 2,
            n => unreachable!("base_to_ctx2: {n}"),
        }
    }

    pub(crate) fn handle_depth(&mut self, bchrom: &str, bpos: i64) {
        if let Some(ref mut depth) = self.depth {
            depth.flush_before(bchrom, bpos);
        }
    }

    pub(crate) fn increment(
        &mut self,
        a: &Record,
        b: &Record,
        a_cigar: &CigarStringView,
        ctx: &ChromCtx<'_>,
    ) {
        let b_cigar = b.cigar();
        let a_qual = a.qual();
        let b_qual = b.qual();
        overlap_pieces_into(
            a_cigar,
            &b_cigar,
            a_qual,
            b_qual,
            true,
            &mut self.pieces_scratch,
        );
        if self.pieces_scratch.is_empty() {
            return;
        }
        let a_seq = a.seq();
        let b_seq = b.seq();
        let bins = self.counts.cnts.dim().2;
        let max_len = a_seq.len().max(b_seq.len());
        if max_len / self.bin_size as usize >= bins {
            panic!(
                "index out of bounds: specify a --max-read-length of at least {}",
                max_len
            );
        }

        let ma = Mate::new(a, a_cigar, bins);
        let mb = Mate::new(b, &b_cigar, bins);
        let min_base_qual = self.min_base_qual;
        let min_map_qual = self.min_map_qual;
        let bin_size = self.bin_size;

        if has_indel(a_cigar) || has_indel(&b_cigar) {
            let indel_errors =
                indel_error_pieces_with_qpos(a_cigar, &b_cigar, a_qual, b_qual, min_base_qual);
            for ie in &indel_errors {
                let c = &ie.coord;
                // include the event if any of it overlaps with the include tree.
                if let Some(t) = ctx.include_tree {
                    if t.count(c.start, c.stop) == 0 {
                        continue;
                    }
                }
                // exclude the event if any of it overlaps with the exclude tree.
                if let Some(t) = ctx.exclude_tree {
                    if t.count(c.start, c.stop) != 0 {
                        continue;
                    }
                }

                let len = match c.indel_type {
                    IndelType::Insertion(len) => len as i32,
                    IndelType::Deletion(len) => -(len as i32),
                    IndelType::NotIndel => 0,
                };

                if let Some(lf) = self.lua_filter.as_ref().filter(|lf| lf.has_base_filter()) {
                    let (rec, seq) = if ie.from_a { (a, &a_seq) } else { (b, &b_seq) };
                    let base_byte = if ie.qpos < seq.len() { seq[ie.qpos] } else { 0 };
                    match lf.skip_base(
                        rec,
                        BasePos {
                            qpos: ie.qpos,
                            gpos: c.start,
                            seq_len: seq.len(),
                            bq: c.qual,
                            base: base_byte,
                            is_reverse: rec.is_reverse(),
                            is_indel: true,
                            indel_length: len,
                        },
                    ) {
                        Ok(true) => continue,
                        Ok(false) => (),
                        Err(e) => {
                            if !self.logged_base_err {
                                log::error!(
                                    "error evaluating base lua expression for indel (further errors suppressed): {e}"
                                );
                                self.logged_base_err = true;
                            }
                        }
                    }
                }

                let p = Position {
                    tid: a.tid() as u16,
                    pos: c.start,
                    bq_bin: Counts::qual_to_bin(c.qual),
                };

                let hps = ctx.hp_tree.as_ref().map(|t| {
                    t.find(
                        c.start.max(MAX_HP_DIST as u32) - MAX_HP_DIST as u32,
                        c.stop + MAX_HP_DIST as u32,
                    )
                    .collect::<Vec<_>>()
                });

                let hp_dist_a =
                    hp::hp_distance(hps.as_deref(), c.start, ma.start, ma.end, ma.strand);
                let hp_dist_b =
                    hp::hp_distance(hps.as_deref(), c.start, mb.start, mb.end, mb.strand);
                let indel_hp_dist =
                    min_abs_hp_distance(hp_dist_a, hp_dist_b).unwrap_or(MAX_HP_DIST + 1);

                *self
                    .counts
                    .indel_error_positions
                    .entry((p, len, indel_hp_dist))
                    .or_insert(0) += 1;
            }
        }

        let mut genome_pos = u32::MAX;
        // Take `pieces_scratch` out so we can borrow other `self` fields freely inside the loop,
        // then put it back to preserve its heap capacity across pairs.
        let pieces = std::mem::take(&mut self.pieces_scratch);
        debug_assert_eq!(
            self.counts.cnts.strides(),
            &[
                (2 * bins * CNTS_BIN_STRIDE) as isize,
                (bins * CNTS_BIN_STRIDE) as isize,
                CNTS_BIN_STRIDE as isize,
                64,
                32,
                1,
            ]
        );
        // Hoist `cnts.as_slice_mut()` outside the loops so ndarray's 6D layout check runs once per pair.
        let cnts_slice = self.counts.cnts.as_slice_mut().unwrap();
        let bin_off_table = &self.bin_off_table;

        for [a_chunk, b_chunk, g_chunk] in pieces.iter() {
            // Query include/exclude trees once per chunk instead of per base; skip the whole
            // chunk (including Lua and homopolymer queries) if no base can pass.
            if self
                .inc_filter
                .prepare(ctx.include_tree, g_chunk.start, g_chunk.stop)
                == ChunkFilter::SkipAll
                || self
                    .exc_filter
                    .prepare(ctx.exclude_tree, g_chunk.start, g_chunk.stop)
                    == ChunkFilter::SkipAll
            {
                continue;
            }

            let chunk_len = (a_chunk.stop - a_chunk.start) as usize;
            let use_mask = if let Some(lf) =
                self.lua_filter.as_ref().filter(|lf| lf.has_base_filter())
            {
                let (a_seq_slice, b_seq_slice) = if lf.needs_base_char() {
                    self.a_seq_scratch.clear();
                    self.b_seq_scratch.clear();
                    for ai in a_chunk.start..a_chunk.stop {
                        self.a_seq_scratch
                            .push(unsafe { a_seq.decoded_base_unchecked(ai as usize) });
                    }
                    for bi in b_chunk.start..b_chunk.stop {
                        self.b_seq_scratch
                            .push(unsafe { b_seq.decoded_base_unchecked(bi as usize) });
                    }
                    (
                        Some(self.a_seq_scratch.as_slice()),
                        Some(self.b_seq_scratch.as_slice()),
                    )
                } else {
                    (None, None)
                };

                match lf.eval_chunk(
                    ChunkEval {
                        a,
                        b,
                        a_q0: a_chunk.start,
                        b_q0: b_chunk.start,
                        g0: g_chunk.start,
                        len: chunk_len,
                        min_bq: min_base_qual,
                        a_quals: &a_qual[a_chunk.start as usize..a_chunk.stop as usize],
                        b_quals: &b_qual[b_chunk.start as usize..b_chunk.stop as usize],
                        a_seq: a_seq_slice,
                        b_seq: b_seq_slice,
                    },
                    &mut self.base_skip_mask,
                ) {
                    Ok(ChunkSkip::All) => continue,
                    Ok(ChunkSkip::None) => false,
                    Ok(ChunkSkip::Mask) => true,
                    Err(e) => {
                        if !self.logged_base_err {
                            log::error!(
                                "error evaluating base lua expression (further errors suppressed): {e}"
                            );
                            self.logged_base_err = true;
                        }
                        false
                    }
                }
            } else {
                false
            };

            // we want to limit to the bounds of the read. since homopolymers outside of the read won't affect it.
            let eps = 1;
            let g_start = (g_chunk.start.max(MAX_HP_DIST as u32) - MAX_HP_DIST as u32)
                .max(ma.start + eps)
                .max(mb.start + eps);
            let g_stop = (g_chunk.stop + MAX_HP_DIST as u32)
                .min(ma.end - eps)
                .min(mb.end - eps);
            let hps: Option<Vec<_>> = if g_start <= g_stop {
                ctx.hp_tree.as_ref().map(|t| t.find(g_start, g_stop).collect())
            } else {
                None
            };

            // Monotonic cursor into `hps`: since `genome_pos` strictly increases across the chunk
            // and `hps` is sorted by start, we never rescan homopolymers > MAX_HP_DIST behind us.
            let mut hp_cursor: usize = 0;

            for (ai, bi) in std::iter::zip(a_chunk.start..a_chunk.stop, b_chunk.start..b_chunk.stop)
            {
                let aq = a_qual[ai as usize];
                if aq < min_base_qual {
                    continue;
                }
                let bq = b_qual[bi as usize];
                if bq < min_base_qual {
                    continue;
                }
                let offset = ai - a_chunk.start;
                if use_mask && self.base_skip_mask[offset as usize] != 0 {
                    continue;
                }
                genome_pos = g_chunk.start + offset;

                if !self.inc_filter.passes(genome_pos) || !self.exc_filter.passes(genome_pos) {
                    continue;
                }

                // Read raw 4-bit BAM base nibbles first so matching bases (>99%) avoid ASCII
                // decoding (`ENC_TO_CTX2` is a 16-byte LUT). Skip N / ambiguous IUPAC bases
                // (ctx2 == 2) before updating depth, counts, or mismatches (`context_lookup`
                // only handles A/C/G/T pairs and `cnts` has only 2 context slots).
                let a_enc = unsafe { a_seq.encoded_base_unchecked(ai as usize) };
                let b_enc = unsafe { b_seq.encoded_base_unchecked(bi as usize) };
                let mut a_ctx2 = ENC_TO_CTX2[a_enc as usize];
                let mut b_ctx2 = ENC_TO_CTX2[b_enc as usize];
                if (a_ctx2 | b_ctx2) >= 2 {
                    continue;
                }

                let aq = Counts::qual_to_bin(aq);
                let bq = Counts::qual_to_bin(bq);

                if let Some(ref mut depth) = self.depth {
                    depth.add(genome_pos, aq, bq);
                }

                let (a_hp_dist, b_hp_dist) = if let Some(ref hps_vec) = hps {
                    while hp_cursor < hps_vec.len()
                        && hps_vec[hp_cursor].stop + (MAX_HP_DIST as u32) < genome_pos
                    {
                        hp_cursor += 1;
                    }
                    let slice = &hps_vec[hp_cursor..];
                    let ah = hp::hp_distance(Some(slice), genome_pos, ma.start, ma.end, ma.strand)
                        .map(|d| (d + MAX_HP_DIST) as usize)
                        .unwrap_or(NO_HP_BIN);
                    let bh = hp::hp_distance(Some(slice), genome_pos, mb.start, mb.end, mb.strand)
                        .map(|d| (d + MAX_HP_DIST) as usize)
                        .unwrap_or(NO_HP_BIN);
                    (ah, bh)
                } else {
                    (NO_HP_BIN, NO_HP_BIN)
                };

                let a_bin_off = bin_off_table[ai as usize];
                let b_bin_off = bin_off_table[bi as usize];

                if a_enc == b_enc {
                    cnts_slice[cnts_flat_idx(ma.row_base, a_bin_off, aq, a_ctx2, a_hp_dist)] += 1;
                    cnts_slice[cnts_flat_idx(mb.row_base, b_bin_off, bq, a_ctx2, b_hp_dist)] += 1;
                    self.counts.matches += 1;
                    continue;
                }

                let a_base = ENC_TO_BASE[a_enc as usize];
                let b_base = ENC_TO_BASE[b_enc as usize];

                let a_index = [
                    ma.read12,
                    ma.rev_idx,
                    (ai / bin_size) as usize,
                    aq as usize,
                    a_ctx2,
                    a_hp_dist,
                ];
                let b_index = [
                    mb.read12,
                    mb.rev_idx,
                    (bi / bin_size) as usize,
                    bq as usize,
                    b_ctx2,
                    b_hp_dist,
                ];

                self.counts.mismatches += 1;

                let real_base = if let Some(ref mut pileup) = self.pileup {
                    let mut base_counts =
                        pileup.query(a.tid(), genome_pos, min_map_qual, min_base_qual);
                    let am = argmax(&base_counts).expect("error selecting maximum index");
                    // check that the 2nd most common base is very low frequency, otherwise might be a het.
                    base_counts.sort();
                    let cmax = base_counts[4];
                    // if 2nd most common base is more than 50% of first, then we don't know which is right.
                    if base_counts[3] as f64 / cmax as f64 > 0.5 {
                        log::debug!(
                            "skipping due to unknown truth given base_counts {:?} at pos:{}:{}",
                            base_counts,
                            ctx.chrom,
                            genome_pos
                        );
                        continue;
                    }
                    ['A', 'C', 'G', 'T', 'N'][am]
                } else {
                    ctx.chrom_seq
                        .as_deref()
                        .expect("chrom_seq must be loaded when pileup is disabled")
                        [genome_pos as usize] as char
                };
                if real_base == 'N' {
                    log::warn!("got 'N' for {}:{genome_pos}. skipping", ctx.chrom);
                    let pos = Position {
                        tid: a.tid() as u16,
                        pos: genome_pos,
                        // we don't know the bq, but assume it's the min. this very rarely happens so doesn't affect results.
                        bq_bin: aq.min(bq),
                    };
                    let context_counts = self.counts.error_positions.entry(pos).or_insert([0; 7]);
                    context_counts[6] += 1;
                    continue;
                }

                let err_index = if a_base == real_base as u8 {
                    // b is the error; correct b's context to the true base before incrementing cnts.
                    let mut index = b_index;
                    b_ctx2 = a_ctx2;
                    index[4] = context_lookup(a_base, b_base);
                    index
                } else if b_base == real_base as u8 {
                    // a is the error; correct a's context to the true base before incrementing cnts.
                    let mut index = a_index;
                    a_ctx2 = b_ctx2;
                    index[4] = context_lookup(b_base, a_base);
                    index
                } else {
                    // can't determine which is error base.
                    let pos = Position {
                        tid: a.tid() as u16,
                        pos: genome_pos,
                        // we don't know the bq, but assume it's the min. this very rarely happens so doesn't affect results.
                        bq_bin: aq.min(bq),
                    };
                    log::debug!(
                        "bases mismatches between reads and neither matches reference at pos:{}:{}. adding N",
                        ctx.chrom,
                        genome_pos
                    );
                    let context_counts = self.counts.error_positions.entry(pos).or_insert([0; 7]);
                    context_counts[6] += 1;
                    continue;
                };

                let pos = Position {
                    tid: a.tid() as u16,
                    pos: genome_pos,
                    bq_bin: err_index[3] as u8,
                };
                let context_idx = err_index[4];
                let context_counts = self.counts.error_positions.entry(pos).or_insert([0; 7]);
                context_counts[context_idx] += 1;

                cnts_slice[cnts_flat_idx(ma.row_base, a_bin_off, aq, a_ctx2, a_hp_dist)] += 1;
                cnts_slice[cnts_flat_idx(mb.row_base, b_bin_off, bq, b_ctx2, b_hp_dist)] += 1;

                self.counts.errs[err_index] += 1;
            }
        }
        self.pieces_scratch = pieces;
        if genome_pos != u32::MAX {
            if let Some(ref mut depth) = self.depth {
                // Flush positions strictly before max(a_start, b_start): in coordinate-sorted BAMs,
                // no subsequent overlapping pair can cover positions earlier than the second mate's start.
                depth.flush_before(ctx.chrom, ma.start.max(mb.start) as i64);
            }
        }
    }
}

lazy_static! {
    pub(crate) static ref CONTEXT_TO_CONTEXT2: [[char; 2]; 7] = [
        ['A', 'C'],
        ['A', 'G'],
        ['A', 'T'],
        ['C', 'A'],
        ['C', 'G'],
        ['C', 'T'],
        ['N', 'N'],
    ];
    pub(crate) static ref Q_LOOKUP: [&'static str; 5] = ["0-5", "05-19", "20-36", "37-59", "60+"];
    pub(crate) static ref REVERSE_Q_LOOKUP: HashMap<&'static str, u8> = HashMap::from([
        ("0-5", 0),
        ("05-19", 1),
        ("20-36", 2),
        ("37-59", 3),
        ("60+", 4),
    ]);
}

pub(crate) fn filter_read(r: &Record) -> bool {
    r.tid() == r.mtid()
        && r.tid() >= 0
        && !r.is_unmapped()
        && !r.is_mate_unmapped()
        && (r.pos() - r.mpos()).abs() < 1000
        && !r.is_supplementary()
        && !r.is_secondary()
        && !r.is_duplicate()
        && !r.is_quality_check_failed()
}

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Coordinates {
    pub start: u32,
    pub stop: u32,
    pub indel_type: IndelType,
    pub qual: u8,
}

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Clone)]
pub enum IndelType {
    Insertion(u32),
    Deletion(u32),
    NotIndel,
}

#[inline(always)]
fn is_insertion(a: Cigar) -> bool {
    matches!(a, Cigar::Ins(_))
}

#[inline(always)]
fn query(a: Cigar) -> i64 {
    match a {
        Cigar::Match(n) | Cigar::SoftClip(n) | Cigar::Ins(n) | Cigar::Diff(n) | Cigar::Equal(n) => {
            n as i64
        }
        _ => 0,
    }
}
#[inline(always)]
fn reference(a: Cigar) -> i64 {
    match a {
        Cigar::Match(n) | Cigar::Del(n) | Cigar::Diff(n) | Cigar::Equal(n) | Cigar::RefSkip(n) => {
            n as i64
        }
        _ => 0,
    }
}

struct IndelError {
    coord: Coordinates,
    qpos: usize,
    from_a: bool,
}

fn indel_coords_with_qpos(
    cig: &CigarStringView,
    genomic_min: u32,
    genomic_max: u32,
    base_quals: &[u8],
) -> Vec<(Coordinates, usize)> {
    let mut result: Vec<(Coordinates, usize)> = Vec::new();
    let mut start: u32 = cig.pos() as u32;
    let mut read_i = 0;

    for c in cig {
        if start > genomic_max {
            break;
        }
        if start + (reference(*c) as u32) < genomic_min {
            start += reference(*c) as u32;
            read_i += query(*c) as usize;
            continue;
        }
        // TODO: handle partial overlap of start with the current cigar op.
        match c {
            Cigar::Ins(l) => {
                result.push((
                    Coordinates {
                        start,
                        stop: start + 1,
                        indel_type: IndelType::Insertion(*l ),
                        qual: base_quals[read_i],
                    },
                    read_i,
                ));
            }
            Cigar::Del(d) => {
                result.push((
                    Coordinates {
                        start,
                        stop: start + *d,
                        indel_type: IndelType::Deletion(*d ),
                        qual: base_quals[read_i],
                    },
                    read_i,
                ));
            }
            _ => {}
        }
        start += reference(*c) as u32;
        read_i += query(*c) as usize;
    }
    result
}

fn find_non_exact_with_qpos(
    a_indel_coords: &[(Coordinates, usize)],
    b_indel_coords: &[(Coordinates, usize)],
    result: &mut Vec<IndelError>,
    min_base_qual: u8,
    from_a: bool,
) {
    for (a, qpos) in a_indel_coords {
        if a.qual <= min_base_qual {
            continue;
        }
        match b_indel_coords.binary_search_by(|(b, _)| b.cmp(a)) {
            Ok(_) => {}
            Err(bi) => {
                // we check base-qual (first base) of b as well.
                // this is a bit weird, but ensures that at least both
                // reads were confident at this site.
                if bi < b_indel_coords.len() && b_indel_coords[bi].0.qual < min_base_qual {
                    continue;
                }
                let bq = if bi < b_indel_coords.len() {
                    b_indel_coords[bi].0.qual
                } else {
                    u8::MAX
                };
                // any non-exact matches are errors
                result.push(IndelError {
                    coord: Coordinates {
                        start: a.start,
                        stop: a.stop,
                        indel_type: a.indel_type.clone(),
                        qual: a.qual.min(bq),
                    },
                    qpos: *qpos,
                    from_a,
                });
            }
        }
    }
}

/// Report genomic coordinates (and query positions) of indels that do not match between the reads.
fn indel_error_pieces_with_qpos(
    a: &CigarStringView,
    b: &CigarStringView,
    a_qual: &[u8],
    b_qual: &[u8],
    min_base_qual: u8,
) -> Vec<IndelError> {
    let aend = a.end_pos() as u32;
    let bend = b.end_pos() as u32;
    if aend <= b.pos() as u32 || bend <= a.pos() as u32 {
        return vec![];
    }
    let a_indel_coords = indel_coords_with_qpos(a, b.pos() as u32, bend, a_qual);
    let b_indel_coords = indel_coords_with_qpos(b, a.pos() as u32, aend, b_qual);
    if a_indel_coords.is_empty() && b_indel_coords.is_empty() {
        return vec![];
    }

    let mut result: Vec<IndelError> = Vec::new();
    find_non_exact_with_qpos(&a_indel_coords, &b_indel_coords, &mut result, min_base_qual, true);
    find_non_exact_with_qpos(&b_indel_coords, &a_indel_coords, &mut result, min_base_qual, false);
    result
}

/// Report genomic coordiantes of bases that do not match between the reads.
#[cfg(test)]
fn indel_error_pieces(
    a: &CigarStringView,
    b: &CigarStringView,
    a_qual: &[u8],
    b_qual: &[u8],
    min_base_qual: u8,
) -> Vec<Coordinates> {
    indel_error_pieces_with_qpos(a, b, a_qual, b_qual, min_base_qual)
        .into_iter()
        .map(|ie| ie.coord)
        .collect()
}

/// Return mapped parts of each read that overlap the other into `result` (reusing its capacity).
/// Each entry holds `[a_read_coords, b_read_coords, genome_coords]`.
fn overlap_pieces_into(
    a: &CigarStringView,
    b: &CigarStringView,
    a_qual: &[u8],
    b_qual: &[u8],
    skip_insertions: bool,
    result: &mut Vec<[Coordinates; 3]>,
) {
    result.clear();
    // Fast path for the common case (~95%+ of short-read pairs): both reads are a single
    // ungapped alignment op (M/=/X), so the overlap is one contiguous interval with no loop.
    if a.len() == 1
        && b.len() == 1
        && matches!(a[0], Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_))
        && matches!(b[0], Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_))
    {
        let a_start = a.pos();
        let b_start = b.pos();
        let a_end = a_start + a[0].len() as i64;
        let b_end = b_start + b[0].len() as i64;
        let g_start = a_start.max(b_start);
        let g_stop = a_end.min(b_end);
        if g_stop <= g_start {
            return;
        }
        let glen = (g_stop - g_start) as u32;
        let a_over = (g_start - a_start) as u32;
        let b_over = (g_start - b_start) as u32;
        result.push([
            Coordinates {
                start: a_over,
                stop: a_over + glen,
                indel_type: IndelType::NotIndel,
                qual: a_qual[a_over as usize],
            },
            Coordinates {
                start: b_over,
                stop: b_over + glen,
                indel_type: IndelType::NotIndel,
                qual: b_qual[b_over as usize],
            },
            Coordinates {
                start: g_start as u32,
                stop: g_stop as u32,
                indel_type: IndelType::NotIndel,
                qual: a_qual[a_over as usize],
            },
        ]);
        return;
    }

    let aend = a.end_pos();
    let bend = b.end_pos();
    if aend <= b.pos() || bend <= a.pos() {
        return;
    }

    let mut ai: usize = 0;
    let mut bi: usize = 0;
    let mut a_genome_pos = a.pos();
    let mut b_genome_pos = b.pos();
    let mut a_read_pos = 0i64;
    let mut b_read_pos = 0i64;
    while ai < a.len() && bi < b.len() {
        let a_genome_stop = a_genome_pos + reference(a[ai.min(a.len() - 1)]);
        let b_genome_stop = b_genome_pos + reference(b[bi.min(b.len() - 1)]);
        if a_genome_stop < b_genome_pos {
            if ai < a.len() {
                a_genome_pos += reference(a[ai]);
                a_read_pos += query(a[ai]);
                ai += 1;
            }
        } else if b_genome_stop < a_genome_pos {
            if bi < b.len() {
                b_genome_pos += reference(b[bi]);
                b_read_pos += query(b[bi]);
                bi += 1;
            }
        } else {
            // we have some overlap.
            // if they both consume query, we can append to our result.
            let aop = a[ai.min(a.len() - 1)];
            let bop = b[bi.min(b.len() - 1)];
            if query(aop) > 0 && query(bop) > 0 {
                let genome_start = a_genome_pos.max(b_genome_pos);
                let genome_stop = a_genome_stop.min(b_genome_stop);

                let mut glen = genome_stop - genome_start;
                if glen == 0 && !skip_insertions {
                    // if they are both the same insertion, then we will evaluate.
                    // otherwise, we can not.
                    if aop == bop && is_insertion(aop) {
                        glen = aop.len() as i64;
                    }
                }

                // if glen is 0, we didn't consume any reference, but can have, e.g. both deletions.
                if glen > 0 {
                    //let a_over = aop.len() as i64 - (genome_start - a_genome_pos);
                    //let b_over = bop.len() as i64 - (genome_start - b_genome_pos);

                    // glen can be 0 if, e.g. both reads end with soft-clip.
                    let a_over = genome_start - a_genome_pos;
                    let b_over = genome_start - b_genome_pos;

                    result.push([
                        Coordinates {
                            start: (a_read_pos + a_over) as u32,
                            stop: (a_read_pos + a_over + glen) as u32,
                            indel_type: match aop {
                                Cigar::Ins(l) => IndelType::Insertion(l),
                                Cigar::Del(l) => IndelType::Deletion(l),
                                _ => IndelType::NotIndel,
                            },
                            qual: a_qual[a_read_pos as usize + a_over as usize],
                        },
                        Coordinates {
                            start: (b_read_pos + b_over) as u32,
                            stop: (b_read_pos + b_over + glen) as u32,
                            indel_type: match bop {
                                Cigar::Ins(l) => IndelType::Insertion(l),
                                Cigar::Del(l) => IndelType::Deletion(l),
                                _ => IndelType::NotIndel,
                            },
                            qual: b_qual[b_read_pos as usize + b_over as usize],
                        },
                        Coordinates {
                            start: genome_start as u32,
                            stop: genome_stop as u32,
                            indel_type: match aop {
                                Cigar::Ins(l) => IndelType::Insertion(l),
                                Cigar::Del(l) => IndelType::Deletion(l),
                                _ => IndelType::NotIndel,
                            },
                            qual: a_qual[a_read_pos as usize + a_over as usize],
                        },
                    ])
                }
            }
            // we had some overlap. now we increment the lowest genome pos by end.
            if a_genome_stop <= b_genome_stop && ai < a.len() {
                a_genome_pos += reference(a[ai]);
                a_read_pos += query(a[ai]);
                ai += 1;
            }
            if b_genome_stop <= a_genome_stop && bi < b.len() {
                b_genome_pos += reference(b[bi]);
                b_read_pos += query(b[bi]);
                bi += 1;
            }
        }
    }
}

#[cfg(test)]
fn overlap_pieces(
    a: &CigarStringView,
    b: &CigarStringView,
    a_qual: &[u8],
    b_qual: &[u8],
    skip_insertions: bool,
) -> Vec<[Coordinates; 3]> {
    let mut result = Vec::new();
    overlap_pieces_into(a, b, a_qual, b_qual, skip_insertions, &mut result);
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{Cigar, CigarString};

    #[test]
    fn test_different_alignments() {
        let a = CigarString(vec![Cigar::Match(5), Cigar::Ins(3), Cigar::Match(5)]).into_view(0);
        let b = CigarString(vec![Cigar::Match(13)]).into_view(0);
        let a_bqs = vec![30u8; 20];
        let b_bqs = vec![30u8; 20];
        let r = overlap_pieces(&a, &b, &a_bqs, &b_bqs, true);
        dbg!(&r);
    }

    #[test]
    fn test_same_insertion() {
        let a = CigarString(vec![Cigar::Match(10), Cigar::Ins(8), Cigar::Match(10)]).into_view(5);
        let b = CigarString(vec![Cigar::Match(10), Cigar::Ins(8), Cigar::Match(10)]).into_view(5);
        let a_bqs = vec![30u8; 20];
        let b_bqs = vec![30u8; 20];
        let r = overlap_pieces(&a, &b, &a_bqs, &b_bqs, false);
        dbg!(&r);
    }

    #[test]
    fn test_many_contained() {
        let a = CigarString(vec![Cigar::Match(100)]).into_view(10);
        let b = CigarString(vec![
            Cigar::Match(10),
            Cigar::Match(11),
            Cigar::Match(12),
            Cigar::Match(13),
        ])
        .into_view(0);
        let a_bqs = vec![30u8; 200];
        let b_bqs = vec![30u8; 200];
        let r = overlap_pieces(&a, &b, &a_bqs, &b_bqs, true);
        let expected = [
            [
                Coordinates {
                    start: 0,
                    stop: 11,
                    indel_type: IndelType::NotIndel,
                    qual: 30,
                },
                Coordinates {
                    start: 10,
                    stop: 21,
                    indel_type: IndelType::NotIndel,
                    qual: 30,
                },
                Coordinates {
                    start: 10,
                    stop: 21,
                    indel_type: IndelType::NotIndel,
                    qual: 30,
                },
            ],
            [
                Coordinates {
                    start: 11,
                    stop: 23,
                    indel_type: IndelType::NotIndel,
                    qual: 30,
                },
                Coordinates {
                    start: 21,
                    stop: 33,
                    indel_type: IndelType::NotIndel,
                    qual: 30,
                },
                Coordinates {
                    start: 21,
                    stop: 33,
                    indel_type: IndelType::NotIndel,
                    qual: 30,
                },
            ],
            [
                Coordinates {
                    start: 23,
                    stop: 36,
                    indel_type: IndelType::NotIndel,
                    qual: 30,
                },
                Coordinates {
                    start: 33,
                    stop: 46,
                    indel_type: IndelType::NotIndel,
                    qual: 30,
                },
                Coordinates {
                    start: 33,
                    stop: 46,
                    indel_type: IndelType::NotIndel,
                    qual: 30,
                },
            ],
        ];
        assert_eq!(r, expected);
    }

    #[test]
    fn test_simple_overlap() {
        let a = CigarString(vec![
            Cigar::Match(10),
            Cigar::Match(80),
            Cigar::SoftClip(10),
        ])
        .into_view(8);
        let b = CigarString(vec![
            Cigar::Match(70),
            Cigar::Match(40),
            Cigar::SoftClip(10),
        ])
        .into_view(5);

        let a_bqs = vec![30u8; 100];
        let b_bqs = vec![30u8; 100];
        let r = overlap_pieces(&a, &b, &a_bqs, &b_bqs, true);

        let expected = [
            [
                Coordinates {
                    start: 0,
                    stop: 10,
                    indel_type: IndelType::NotIndel,
                    qual: 30,
                },
                Coordinates {
                    start: 3,
                    stop: 13,
                    indel_type: IndelType::NotIndel,
                    qual: 30,
                },
                Coordinates {
                    start: 8,
                    stop: 18,
                    indel_type: IndelType::NotIndel,
                    qual: 30,
                },
            ],
            [
                Coordinates {
                    start: 10,
                    stop: 67,
                    indel_type: IndelType::NotIndel,
                    qual: 30,
                },
                Coordinates {
                    start: 13,
                    stop: 70,
                    indel_type: IndelType::NotIndel,
                    qual: 30,
                },
                Coordinates {
                    start: 18,
                    stop: 75,
                    indel_type: IndelType::NotIndel,
                    qual: 30,
                },
            ],
            [
                Coordinates {
                    start: 67,
                    stop: 90,
                    indel_type: IndelType::NotIndel,
                    qual: 30,
                },
                Coordinates {
                    start: 70,
                    stop: 93,
                    indel_type: IndelType::NotIndel,
                    qual: 30,
                },
                Coordinates {
                    start: 75,
                    stop: 98,
                    indel_type: IndelType::NotIndel,
                    qual: 30,
                },
            ],
        ];

        assert_eq!(r, expected);
    }

    #[test]
    fn test_size() {
        assert_eq!(std::mem::size_of::<Position>(), 8);
    }

    #[test]
    fn test_min_abs_hp_distance() {
        // Both distances available - choose minimum absolute value
        assert_eq!(min_abs_hp_distance(Some(-1), Some(2)), Some(-1));
        assert_eq!(min_abs_hp_distance(Some(2), Some(-1)), Some(-1));
        assert_eq!(min_abs_hp_distance(Some(-3), Some(-2)), Some(-2));
        assert_eq!(min_abs_hp_distance(Some(3), Some(2)), Some(2));

        // Equal absolute values - choose first one
        assert_eq!(min_abs_hp_distance(Some(-2), Some(2)), Some(-2));
        assert_eq!(min_abs_hp_distance(Some(2), Some(-2)), Some(2));

        // Only one distance available
        assert_eq!(min_abs_hp_distance(Some(5), None), Some(5));
        assert_eq!(min_abs_hp_distance(None, Some(-3)), Some(-3));

        // No distances available
        assert_eq!(min_abs_hp_distance(None, None), None);
    }

    #[test]
    fn test_indel_error_pieces() {
        let cigar_a =
            CigarString(vec![Cigar::Match(10), Cigar::Del(2), Cigar::Match(5)]).into_view(5);
        let cigar_b =
            CigarString(vec![Cigar::Match(10), Cigar::Del(3), Cigar::Match(4)]).into_view(10);
        let expected = vec![
            Coordinates {
                start: 15,
                stop: 17,
                indel_type: IndelType::Deletion(2),
                qual: 30,
            },
            Coordinates {
                start: 20,
                stop: 23,
                indel_type: IndelType::Deletion(3),
                qual: 30,
            },
        ];
        let a_bqs = vec![30u8; 20];
        let b_bqs = vec![30u8; 20];
        assert_eq!(
            indel_error_pieces(&cigar_a, &cigar_b, &a_bqs, &b_bqs, 15),
            expected
        );
    }

    #[test]
    fn test_indel_error_pieces_overlap() {
        let cigar_a =
            CigarString(vec![Cigar::Match(10), Cigar::Del(3), Cigar::Match(5)]).into_view(10);
        let cigar_b =
            CigarString(vec![Cigar::Match(10), Cigar::Ins(3), Cigar::Match(5)]).into_view(10);
        let expected = vec![
            Coordinates {
                start: 20,
                stop: 23,
                indel_type: IndelType::Deletion(3),
                qual: 30,
            },
            Coordinates {
                start: 20,
                stop: 21,
                indel_type: IndelType::Insertion(3),
                qual: 30,
            },
        ];
        let a_bqs = vec![30u8; 20];
        let b_bqs = vec![30u8; 20];
        assert_eq!(
            indel_error_pieces(&cigar_a, &cigar_b, &a_bqs, &b_bqs, 10),
            expected
        );
    }

    #[test]
    fn test_indel_error_quals() {
        let cigar_a =
            CigarString(vec![Cigar::Match(10), Cigar::Del(3), Cigar::Match(5)]).into_view(10);
        let cigar_b =
            CigarString(vec![Cigar::Match(10), Cigar::Ins(3), Cigar::Match(5)]).into_view(10);
        let expected = vec![];
        let a_bqs = vec![10u8; 20];
        let b_bqs = vec![10u8; 20];
        assert_eq!(
            indel_error_pieces(&cigar_a, &cigar_b, &a_bqs, &b_bqs, 60),
            expected
        );
    }

    #[test]
    fn test_n_bases_skipped() {
        // Test both matching N/N (including qual 65 -> bin 4) and mismatching N/A against truth A:
        // neither should panic, inflate counts/mismatches, or emit denominator depth.
        let mut a = Record::new();
        a.set(
            b"pair1",
            Some(&CigarString(vec![Cigar::Match(5)])),
            b"ANNNA",
            &[35, 35, 65, 35, 35],
        );
        a.set_tid(0);
        a.set_pos(100);
        a.set_flags(0x1 | 0x40); // paired, read1

        let mut b = Record::new();
        b.set(
            b"pair1",
            Some(&CigarString(vec![Cigar::Match(5)])),
            b"ANNAA",
            &[35, 35, 65, 35, 35],
        );
        b.set_tid(0);
        b.set_pos(100);
        b.set_flags(0x1 | 0x80 | 0x10); // paired, read2, reverse

        let mut ref_seq = vec![b'A'; 120];
        ref_seq[100..105].copy_from_slice(b"AAAAA");

        let tmp = std::env::temp_dir().join(format!(
            "fraguracy-n-depth-test-{}.bed.gz",
            std::process::id()
        ));
        let tmp_str = tmp.to_string_lossy().into_owned();

        let mut counts = Counts::new(None, 10, 10, 20, 20, None);
        counts.set_depth_writer(&tmp_str).unwrap();
        let a_cigar = a.cigar();
        let ctx = ChromCtx {
            chrom: "chr1",
            chrom_seq: Some(ref_seq),
            include_tree: None,
            exclude_tree: None,
            hp_tree: None,
        };
        counts.increment(&a, &b, &a_cigar, &ctx);
        counts.handle_depth("chr1", i64::MAX);
        // Drop `counts` so the BGZF writer flushes and closes cleanly before reading back.
        let matches = counts.counts.matches;
        let mismatches = counts.counts.mismatches;
        let cnts_sum = counts.counts.cnts.sum();
        let err_pos_len = counts.counts.error_positions.len();
        drop(counts);

        // Only positions 100 (A/A) and 104 (A/A) should be counted; 101-102 (N/N) and 103 (N/A) skipped.
        assert_eq!(matches, 2);
        assert_eq!(mismatches, 0);
        assert_eq!(cnts_sum, 4);
        assert_eq!(err_pos_len, 0);

        use std::io::Read;
        let mut gz = flate2::read::MultiGzDecoder::new(std::fs::File::open(&tmp).unwrap());
        let mut depth_bed = String::new();
        gz.read_to_string(&mut depth_bed).unwrap();
        let _ = std::fs::remove_file(&tmp);
        let rows: Vec<&str> = depth_bed.lines().collect();
        assert_eq!(
            rows,
            vec![
                "#chrom\tstart\tend\tread1_bq_bin\tread2_bq_bin\tpair-ovl-depth",
                "chr1\t100\t101\t20-36\t20-36\t1",
                "chr1\t104\t105\t20-36\t20-36\t1",
            ]
        );
    }

    #[test]
    fn test_chunk_include_exclude_filters() {
        use rust_lapper::Interval;
        let iv = |start, stop| Interval {
            start,
            stop,
            val: 0u8,
        };

        // Unmerged overlapping & bookended intervals in include_tree, plus partial exclude_tree.
        let inc = Lapper::new(vec![iv(100, 104), iv(103, 106), iv(108, 110)]);
        let exc = Lapper::new(vec![iv(101, 103), iv(102, 104)]);
        // Positions in [100, 110):
        //   included: 100, 101, 102, 103, 104, 105, 108, 109 (8 bases)
        //   excluded: 101, 102, 103
        //   kept:     100, 104, 105, 108, 109 (5 bases)
        let mut a = Record::new();
        a.set(
            b"pair1",
            Some(&CigarString(vec![Cigar::Match(10)])),
            b"ACGTACGTAC",
            &[35; 10],
        );
        a.set_tid(0);
        a.set_pos(100);
        a.set_flags(0x1 | 0x40);

        let mut b = Record::new();
        b.set(
            b"pair1",
            Some(&CigarString(vec![Cigar::Match(10)])),
            b"ACGTACGTAC",
            &[35; 10],
        );
        b.set_tid(0);
        b.set_pos(100);
        b.set_flags(0x1 | 0x80 | 0x10);

        let mut counts = Counts::new(None, 10, 10, 20, 20, None);
        let a_cigar = a.cigar();
        let ctx = ChromCtx {
            chrom: "chr1",
            chrom_seq: None,
            include_tree: Some(&inc),
            exclude_tree: Some(&exc),
            hp_tree: None,
        };
        counts.increment(&a, &b, &a_cigar, &ctx);
        assert_eq!(counts.counts.matches, 5);
        assert_eq!(counts.counts.cnts.sum(), 10);

        // Full-coverage include (via bookended merge) & full-coverage exclude fast paths.
        let mut inc_f = RegionFilter::new(true);
        let mut exc_f = RegionFilter::new(false);
        let full_inc = Lapper::new(vec![iv(100, 105), iv(105, 110)]);
        assert_eq!(inc_f.prepare(Some(&full_inc), 100, 110), ChunkFilter::PassAll);
        assert_eq!(inc_f.prepare(Some(&full_inc), 200, 210), ChunkFilter::SkipAll);
        assert_eq!(exc_f.prepare(Some(&full_inc), 100, 110), ChunkFilter::SkipAll);
        assert_eq!(exc_f.prepare(Some(&full_inc), 200, 210), ChunkFilter::PassAll);
    }
}
