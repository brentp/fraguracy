use bpci::*;
use ndarray::prelude::Array;
use ndarray::Array6;
use rust_htslib::{
    bam::{
        record::{Cigar, CigarStringView},
        IndexedReader, Read, Record,
    },
    bgzf::CompressionLevel,
};
use rust_htslib::{bgzf, faidx};
use rust_lapper::Lapper;
use std::collections::BTreeMap;

use crate::homopolymer as hp;
use std::collections::HashMap;
use std::fmt;
use std::hash::Hash;
use std::io::Write;
use std::rc::Rc;
use std::str;

#[derive(Eq, Hash, PartialEq, Ord, PartialOrd)]
pub(crate) struct Position {
    pub tid: u16,
    pub pos: u32,
    pub bq_bin: u8,
}

/// DepthMap is for a given genome position, the depth at each (aq, bq) pair.
type DepthMap = HashMap<(u8, u8), u32>;
type Length = i32;

pub(crate) const MAX_HP_DIST: i16 = 15;

pub(crate) struct Counts {
    pub(crate) ibam: Option<IndexedReader>,
    //  read, f/r pos, bq, bp, ctx{6} */
    pub(crate) counts: InnerCounts,
    pub(crate) depth: BTreeMap<u32, DepthMap>,
    pub(crate) last_depth_entry: Option<(String, u32, u32, String, String, u32)>,
    pub(crate) depth_writer: Option<bgzf::Writer>,
}

pub(crate) struct InnerCounts {
    // genome_pos
    pub(crate) errs: Array6<u64>,
    //  read, f/r, pos, bq, ctx{2}, hp_dist */
    pub(crate) cnts: Array6<u64>,
    pub(crate) mismatches: u64,
    pub(crate) matches: u64,

    // position -> error count. nice to find sites that are error-prone.
    pub(crate) error_positions: HashMap<Position, [u32; 7]>,
    // position -> indel error counts
    pub(crate) indel_error_positions: HashMap<(Position, Length), u32>,
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
            error_positions: HashMap::new(),
            indel_error_positions: HashMap::new(),
        }
    }
}

impl Counts {
    pub(crate) fn new(ir: Option<IndexedReader>, bins: usize) -> Self {
        Counts {
            ibam: ir,
            counts: InnerCounts::new(bins),
            depth: BTreeMap::new(),
            last_depth_entry: None,
            depth_writer: None,
        }
    }

    pub(crate) fn set_depth_writer(&mut self, path: &str) -> std::io::Result<()> {
        let mut w = bgzf::Writer::from_path_with_level(path, CompressionLevel::Level(1))
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
        w.write_all(b"#chrom\tstart\tend\tread1_bq_bin\tread2_bq_bin\tpair-ovl-depth\n")
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
        self.depth_writer = Some(w);
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
        // this function clears out the BTreeMap of depth entries that are before the current position.
        // it is not called if --no-denominator is specified.
        // it writes out the depth entries as they are popped out of the BTreeMap.
        loop {
            let pos = *(self
                .depth
                .first_key_value()
                .unwrap_or((&u32::MAX, &DepthMap::new()))
                .0);
            if pos == u32::MAX {
                break;
            }
            if (pos as i64) < bpos {
                let depthmap = self.depth.remove(&pos).unwrap();

                for ((aq, bq), dp) in depthmap.iter() {
                    let a_bin = Q_LOOKUP[*aq as usize];
                    let b_bin = Q_LOOKUP[*bq as usize];

                    match &mut self.last_depth_entry {
                        Some((
                            last_chrom,
                            start_pos,
                            last_pos,
                            last_a_bin,
                            last_b_bin,
                            last_dp,
                        )) => {
                            if bchrom == last_chrom
                                && a_bin == last_a_bin
                                && b_bin == last_b_bin
                                && dp == last_dp
                                && *last_pos + 1 == pos
                            {
                                *last_pos = pos;
                            } else {
                                if let Some(writer) = &mut self.depth_writer {
                                    writeln!(
                                        writer,
                                        "{}\t{}\t{}\t{}\t{}\t{}",
                                        last_chrom,
                                        start_pos,
                                        *last_pos + 1,
                                        last_a_bin,
                                        last_b_bin,
                                        last_dp
                                    )
                                    .expect("error writing to bgzf file");
                                }

                                self.last_depth_entry = Some((
                                    bchrom.to_string(),
                                    pos,
                                    pos,
                                    a_bin.to_string(),
                                    b_bin.to_string(),
                                    *dp,
                                ));
                            }
                        }
                        None => {
                            self.last_depth_entry = Some((
                                bchrom.to_string(),
                                pos,
                                pos,
                                a_bin.to_string(),
                                b_bin.to_string(),
                                *dp,
                            ));
                        }
                    }
                }
            } else {
                break;
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn increment<N: AsRef<str> + std::fmt::Debug>(
        &mut self,
        a: Rc<Record>,
        b: Rc<Record>,
        min_base_qual: u8,
        min_map_qual: u8,
        bin_size: u32,
        fasta: &Option<faidx::Reader>,
        chrom: N,
        include_tree: &Option<&Lapper<u32, u8>>,
        exclude_tree: &Option<&Lapper<u32, u8>>,
        hp_tree: &Option<Lapper<u32, u8>>,
    ) {
        let pieces = overlap_pieces(&a.cigar(), &b.cigar(), a.qual(), b.qual(), true);
        if pieces.is_empty() {
            return;
        }
        let a_seq = a.seq();
        let b_seq = b.seq();
        if a_seq.len() / bin_size as usize >= self.counts.cnts.dim().2 {
            panic!(
                "index out of bounds: specify a --max-read-length of at least {}",
                a_seq.len()
            )
        }

        if b_seq.len() / bin_size as usize >= self.counts.cnts.dim().2 {
            panic!(
                "index out of bounds: specify a --max-read-length of at least {}",
                b_seq.len()
            )
        }

        let a_qual = a.qual();
        let b_qual = b.qual();

        let indel_errors =
            indel_error_pieces(&a.cigar(), &b.cigar(), a_qual, b_qual, min_base_qual);
        indel_errors.iter().for_each(|c: &Coordinates| {
            // include the event if any of it overlaps with the include tree.
            if let Some(t) = include_tree {
                if t.count(c.start, c.stop) == 0 {
                    return;
                }
            }
            // exclude the event if any of it overlaps with the exclude tree.
            if let Some(t) = exclude_tree {
                if t.count(c.start, c.stop) != 0 {
                    return;
                }
            }

            let len = match c.indel_type {
                IndelType::Insertion(len) => len as i32,
                IndelType::Deletion(len) => -(len as i32),
                IndelType::NotIndel => 0,
            };

            let p = Position {
                tid: a.tid() as u16,
                pos: c.start,
                bq_bin: Counts::qual_to_bin(c.qual),
            };
            *self
                .counts
                .indel_error_positions
                .entry((p, len))
                .or_insert(0) += 1;
        });

        let mut genome_pos = u32::MAX;
        for [a_chunk, b_chunk, g_chunk] in pieces {
            // we want to limit to the bounds of the read. since homopolymers outside of the read won't affect it.
            let eps = 1;
            let g_start = (g_chunk.start.max(MAX_HP_DIST as u32) - MAX_HP_DIST as u32)
                .max(a.pos() as u32 + eps)
                .max(b.pos() as u32 + eps);
            let g_stop = (g_chunk.stop + MAX_HP_DIST as u32)
                .min(a.cigar().end_pos() as u32 - eps)
                .min(b.cigar().end_pos() as u32 - eps);
            let hps: Option<Vec<_>> = if g_start <= g_stop {
                hp_tree.as_ref().map(|t| t.find(g_start, g_stop).collect())
            } else {
                None
            };

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
                genome_pos = g_chunk.start + (ai - a_chunk.start);

                if let Some(t) = include_tree {
                    if t.count(genome_pos, genome_pos + 1) == 0 {
                        continue;
                    }
                }
                if let Some(t) = exclude_tree {
                    if t.count(genome_pos, genome_pos + 1) != 0 {
                        continue;
                    }
                }

                let aq = Counts::qual_to_bin(aq);
                let bq = Counts::qual_to_bin(bq);

                if self.depth_writer.is_some() {
                    self.depth
                        .entry(genome_pos)
                        .or_default()
                        .entry((aq, bq))
                        .and_modify(|v| *v += 1)
                        .or_insert(1);
                }

                let a_base = unsafe { a_seq.decoded_base_unchecked(ai as usize) };
                let b_base = unsafe { b_seq.decoded_base_unchecked(bi as usize) };

                let a_bin = (ai / bin_size) as usize;
                let b_bin = (bi / bin_size) as usize;

                let a_hp_dist = hp::hp_distance(
                    hps.as_deref(),
                    genome_pos,
                    a.pos() as u32,
                    a.cigar().end_pos() as u32,
                    if a.is_reverse() { -1 } else { 1 },
                )
                .map(|d| (d + MAX_HP_DIST) as usize)
                .unwrap_or((2 * MAX_HP_DIST + 1) as usize);

                let b_hp_dist = hp::hp_distance(
                    hps.as_deref(),
                    genome_pos,
                    b.pos() as u32,
                    b.cigar().end_pos() as u32,
                    if b.is_reverse() { -1 } else { 1 },
                )
                .map(|d| (d + MAX_HP_DIST) as usize)
                .unwrap_or((2 * MAX_HP_DIST + 1) as usize);

                /* read1/2, F/R, pos, mq, bq, ctx, hp_dist */
                let mut a_index = [
                    1 - a.is_first_in_template() as usize, // 0 r1
                    (a.is_reverse() as usize),             //
                    a_bin,
                    aq as usize,
                    // NOTE that this could be an error so we might change this later if we learn a_base is an error
                    Counts::base_to_ctx2(a_base),
                    a_hp_dist,
                ];

                let mut b_index = [
                    1 - b.is_first_in_template() as usize,
                    (b.is_reverse() as usize),
                    b_bin,
                    bq as usize,
                    // NOTE that this could be an error so we might change this later if we learn b_base is an error
                    Counts::base_to_ctx2(b_base),
                    b_hp_dist,
                ];

                if a_base == b_base {
                    // fast path to increment separately here because we must do some extra stuff to error base before incrementing count
                    // if there is an error.
                    self.counts.cnts[a_index] += 1;
                    self.counts.cnts[b_index] += 1;
                    self.counts.matches += 1;
                    continue;
                }

                self.counts.mismatches += 1;
                let mut err = ['X', 'X'];

                let real_base = if self.ibam.is_some() {
                    let mut base_counts = pile(
                        self.ibam.as_mut().unwrap(),
                        a.tid(),
                        genome_pos,
                        min_map_qual,
                        min_base_qual,
                    );
                    let am = argmax(&base_counts).expect("error selecting maximum index");
                    // check that the 2nd most common base is very low frequency, otherwise might be a het.
                    base_counts.sort();
                    let cmax = base_counts[4];
                    // if 3nd most common base is more than 50% of first, then we don't know which is right.
                    if base_counts[3] as f64 / cmax as f64 > 0.5 {
                        log::debug!(
                            "skipping due to unknown truth given base_counts {:?} at pos:{}:{}",
                            base_counts,
                            chrom.as_ref(),
                            genome_pos
                        );
                        continue;
                    }
                    ['A', 'C', 'G', 'T', 'N'][am]
                } else {
                    fasta
                        .as_ref()
                        .unwrap()
                        .fetch_seq(&chrom, genome_pos as usize, genome_pos as usize)
                        .expect("error extracting base")[0] as char
                };
                if real_base == 'N' {
                    log::warn!("got 'N' for {chrom:?}:{genome_pos} base skipping");
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
                    // b is the error
                    let mut index = b_index;
                    b_index[4] = a_index[4]; // we correct this because we want to track the true base
                    index[4] = CONTEXT_LOOKUP[&(a_base, b_base)];
                    err[0] = a_base as char;
                    err[1] = b_base as char;
                    index
                } else if b_base == real_base as u8 {
                    // a is the error
                    let mut index = a_index;
                    a_index[4] = b_index[4]; // we correct this because we want to track the true base
                    index[4] = CONTEXT_LOOKUP[&(b_base, a_base)];
                    err[0] = b_base as char;
                    err[1] = a_base as char;
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
                        chrom.as_ref(),
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

                self.counts.cnts[a_index] += 1;
                self.counts.cnts[b_index] += 1;

                self.counts.errs[err_index] += 1;
                if log::log_enabled!(log::Level::Debug)
                    && unsafe { str::from_utf8_unchecked(a.qname()) }
                        == "A00744:46:HV3C3DSXX:2:2611:30798:35258"
                {
                    log::debug!(
                        "gpos:{}, err:{}->{}, err-index:{:?}, ai:{}, bi:{}, {:?}",
                        genome_pos,
                        /* base_counts, */
                        err[0],
                        err[1],
                        err_index,
                        ai,
                        bi,
                        unsafe { str::from_utf8_unchecked(a.qname()) },
                    );
                }
            }
        }
        if self.depth_writer.is_some() {
            self.handle_depth(chrom.as_ref(), genome_pos as i64);
        }
    }
}

fn pile(
    ibam: &mut IndexedReader,
    tid: i32,
    genome_pos: u32,
    min_map_qual: u8,
    min_base_qual: u8,
) -> [u32; 5] {
    let mut base_counts: [u32; 5] = [0; 5];

    ibam.fetch((tid, genome_pos, genome_pos + 1))
        .expect("Error seeking to genomic position");

    let mut p = ibam.pileup();
    p.set_max_depth(100_000);
    p.filter(|col| col.as_ref().unwrap().pos() == genome_pos)
        .for_each(|col| {
            let col = col.unwrap();

            col.alignments().for_each(|aln| {
                if let Some(qpos) = aln.qpos() {
                    let record = aln.record();
                    // here we want a accurate count, so we skip stuff at either
                    // end of a read (within 3 bases of end)
                    // along with low base-quality and low mapping-quality
                    if qpos < 3 || qpos > record.qual().len() - 4 {
                        return;
                    }
                    if record.mapq() < min_map_qual {
                        return;
                    }
                    if record.qual()[qpos] < min_base_qual {
                        return;
                    }
                    let base_idx = match record.seq()[qpos] as char {
                        'A' => 0,
                        'C' => 1,
                        'G' => 2,
                        'T' => 3,
                        _ => 4,
                    };
                    base_counts[base_idx] += 1;
                }
            });
        });
    base_counts
}

lazy_static! {
    pub(crate) static ref CONTEXT_LOOKUP: HashMap<(u8, u8), usize> = HashMap::from([
        ((b'T', b'G'), 0usize),
        ((b'A', b'C'), 0usize),
        ((b'T', b'C'), 1usize),
        ((b'A', b'G'), 1usize),
        ((b'T', b'A'), 2usize),
        ((b'A', b'T'), 2usize),
        ((b'C', b'A'), 3usize),
        ((b'G', b'T'), 3usize),
        ((b'C', b'G'), 4usize),
        ((b'G', b'C'), 4usize),
        ((b'C', b'T'), 5usize),
        ((b'G', b'A'), 5usize),
        ((b'N', b'N'), 6usize),
    ]);
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

pub(crate) fn filter_read(r: &Rc<Record>) -> bool {
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

fn indel_coords(
    cig: &CigarStringView,
    genomic_min: u32,
    genomic_max: u32,
    base_quals: &[u8],
) -> Vec<Coordinates> {
    let mut result: Vec<Coordinates> = Vec::new();
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
                result.push(Coordinates {
                    start,
                    stop: start + 1,
                    indel_type: IndelType::Insertion(*l as u32),
                    qual: base_quals[read_i],
                });
            }
            Cigar::Del(d) => {
                result.push(Coordinates {
                    start,
                    stop: start + *d,
                    indel_type: IndelType::Deletion(*d as u32),
                    qual: base_quals[read_i],
                });
            }
            _ => {}
        }
        start += reference(*c) as u32;
        read_i += query(*c) as usize;
    }
    result
}

fn find_non_exact(
    a_indel_coords: &[Coordinates],
    b_indel_coords: &[Coordinates],
    result: &mut Vec<Coordinates>,
    min_base_qual: u8,
) {
    for a in a_indel_coords {
        if a.qual <= min_base_qual {
            continue;
        }
        match b_indel_coords.binary_search_by(|b| b.cmp(a)) {
            Ok(_) => {}
            Err(bi) => {
                // we check base-qual (first base) of b as well.
                // this is a bit weird, but ensures that at least both
                // reads were confident at this site.
                if bi < b_indel_coords.len() && b_indel_coords[bi].qual < min_base_qual {
                    continue;
                }
                let bq = if bi < b_indel_coords.len() {
                    b_indel_coords[bi].qual
                } else {
                    u8::MAX
                };
                // any non-exact matches are errors
                result.push(Coordinates {
                    start: a.start,
                    stop: a.stop,
                    indel_type: a.indel_type.clone(),
                    qual: a.qual.min(bq),
                });
            }
        }
    }
}

/// Report genomic coordiantes of bases that do not match between the reads.
fn indel_error_pieces(
    a: &CigarStringView,
    b: &CigarStringView,
    a_qual: &[u8],
    b_qual: &[u8],
    min_base_qual: u8,
) -> Vec<Coordinates> {
    let aend = a.end_pos() as u32;
    let bend = b.end_pos() as u32;
    //let astart = a.pos() + a.leading_softclips();
    //let bstart = b.pos() + b.leading_softclips();
    if aend <= b.pos() as u32 || bend <= a.pos() as u32 {
        return vec![];
    }
    let a_indel_coords = indel_coords(a, b.pos() as u32, bend, a_qual);
    let b_indel_coords = indel_coords(b, a.pos() as u32, aend, b_qual);
    if a_indel_coords.is_empty() && b_indel_coords.is_empty() {
        return vec![];
    }

    let mut result: Vec<Coordinates> = Vec::new();
    find_non_exact(&a_indel_coords, &b_indel_coords, &mut result, min_base_qual);
    find_non_exact(&b_indel_coords, &a_indel_coords, &mut result, min_base_qual);
    result
}

/// Return mapped parts of each read that overlap the other.
/// Returns A, B, genome coordiantes.
fn overlap_pieces(
    a: &CigarStringView,
    b: &CigarStringView,
    a_qual: &[u8],
    b_qual: &[u8],
    skip_insertions: bool,
) -> Vec<[Coordinates; 3]> {
    let aend = a.end_pos();
    let bend = b.end_pos();
    //let astart = a.pos() + a.leading_softclips();
    //let bstart = b.pos() + b.leading_softclips();
    if aend <= b.pos() || bend <= a.pos() {
        return vec![];
    }

    let mut result: Vec<[Coordinates; 3]> = Vec::new();
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
                                Cigar::Ins(l) => IndelType::Insertion(l as u32),
                                Cigar::Del(l) => IndelType::Deletion(l as u32),
                                _ => IndelType::NotIndel,
                            },
                            qual: a_qual[a_read_pos as usize + a_over as usize],
                        },
                        Coordinates {
                            start: (b_read_pos + b_over) as u32,
                            stop: (b_read_pos + b_over + glen) as u32,
                            indel_type: match bop {
                                Cigar::Ins(l) => IndelType::Insertion(l as u32),
                                Cigar::Del(l) => IndelType::Deletion(l as u32),
                                _ => IndelType::NotIndel,
                            },
                            qual: b_qual[b_read_pos as usize + b_over as usize],
                        },
                        Coordinates {
                            start: genome_start as u32,
                            stop: genome_stop as u32,
                            indel_type: match aop {
                                Cigar::Ins(l) => IndelType::Insertion(l as u32),
                                Cigar::Del(l) => IndelType::Deletion(l as u32),
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
}
