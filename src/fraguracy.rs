use ndarray::prelude::{Array, ArrayBase, ArrayView6};
use ndarray::{Array6, ArrayViewMut6};
use rust_htslib::bam::Read;
use rust_htslib::bam::{
    record::{Cigar, CigarStringView},
    IndexedReader, Record,
};
use std::collections::HashMap;
use std::rc::Rc;
use std::str;

pub(crate) struct Counts {
    //  read, pos, mq, bp, ctx{6} */
    pub(crate) ibam: IndexedReader,
    pub(crate) muts: Array6<u64>,
    //  read, pos, mq, bp, ctx{2} */
    pub(crate) cnts: Array6<u64>,
    pub(crate) mismatches: u64,
    pub(crate) matches: u64,
}

fn argmax<T: Ord>(slice: &[T]) -> Option<usize> {
    (0..slice.len()).max_by_key(|i| &slice[*i])
}

impl Counts {
    pub(crate) fn new(ir: IndexedReader) -> Self {
        Counts {
            /*                         read1/2, F/R, pos, mq, bq, ctx */
            ibam: ir,
            cnts: Array::zeros((2, 2, 50, 5, 5, 2)),
            muts: Array::zeros((2, 2, 50, 5, 5, 6)),
            mismatches: 0,
            matches: 0,
        }
    }
    #[inline(always)]
    fn qual_to_bin(q: u8) -> u8 {
        match q {
            0..=5 => 0,
            6..=19 => 1,
            20..=39 => 2,
            40..=59 => 3,
            _ => 4,
        }
    }

    #[inline(always)]
    fn base_to_ctx2(b: u8) -> usize {
        match b as char {
            'A' | 'T' => 0,
            _ => 1,
        }
    }

    pub(crate) fn increment(
        &mut self,
        a: Rc<Record>,
        b: Rc<Record>,
        min_base_qual: u8,
        min_map_qual: u8,
    ) {
        let pieces = overlap_pieces(a.cigar(), b.cigar(), false);
        if pieces.len() == 0 {
            return;
        }
        let a_seq = a.seq();
        let b_seq = b.seq();
        let a_qual = a.qual();
        let b_qual = b.qual();
        let amq = Counts::qual_to_bin(a.mapq());
        let bmq = Counts::qual_to_bin(b.mapq());

        for [a_chunk, b_chunk, g_chunk] in pieces {
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
                let bq = Counts::qual_to_bin(bq);
                let aq = Counts::qual_to_bin(aq);

                let a_base = unsafe { a_seq.decoded_base_unchecked(ai as usize) };
                let b_base = unsafe { b_seq.decoded_base_unchecked(bi as usize) };

                let a_pos = (ai / 3) as usize;
                let b_pos = (bi / 3) as usize;

                /*                         read1/2, F/R, pos, mq, bq, ctx */

                let a_index = [
                    a.is_first_in_template() as usize,
                    1 - (a.is_reverse() as usize),
                    a_pos,
                    amq as usize,
                    aq as usize,
                    Counts::base_to_ctx2(a_base),
                ];

                let b_index = [
                    b.is_first_in_template() as usize,
                    1 - (b.is_reverse() as usize),
                    b_pos,
                    bmq as usize,
                    bq as usize,
                    Counts::base_to_ctx2(b_base),
                ];

                self.cnts[a_index] += 1;
                self.cnts[b_index] += 1;
                self.matches += 1;

                if a_base != b_base {
                    // TODO: pileup and vote to determine error.
                    let genome_pos = g_chunk.start + (ai - a_chunk.start);
                    self.mismatches += 1;
                    let base_counts = pile(
                        &mut self.ibam,
                        a.tid(),
                        genome_pos,
                        min_map_qual,
                        min_base_qual,
                    );
                    let am = argmax(&base_counts).expect("error selecting maximum index");
                    let max_base = ['A', 'C', 'G', 'T'][am];
                    // TODO: check that the 2nd most common base is very low frequency, otherwise might be a het.
                    let mut err = ['A', 'C'];

                    let err_index = if a_base == max_base as u8 {
                        // b is the error
                        let mut index = b_index;
                        index[5] = CONTEXT_LOOKUP[&(a_base, b_base)];
                        err[0] = a_base as char;
                        err[1] = b_base as char;
                        index
                    } else {
                        // a is the error
                        let mut index = a_index;
                        index[5] = CONTEXT_LOOKUP[&(b_base, a_base)];
                        err[0] = b_base as char;
                        err[1] = a_base as char;
                        index
                    };

                    self.muts[err_index] += 1;

                    log::debug!(
                        "gpos: {}, mm: {}, base counts: ACGTN:{:?}, err:{}->{}, err-index:{:?}, ai: {}, bi: {}, {:?}",
                        genome_pos,
                        self.mismatches,
                        base_counts,
                        err[0], err[1],
                        err_index,
                        ai,
                        bi,
                        unsafe { str::from_utf8_unchecked(a.qname()) },
                    );
                }
            }
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
                    let base_idx = match aln.record().seq()[qpos] as char {
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
    return base_counts;
}

lazy_static! {
    pub(crate) static ref CONTEXT_LOOKUP: HashMap<(u8, u8), usize> = HashMap::from([
        (('C' as u8, 'A' as u8), 0usize),
        (('G' as u8, 'T' as u8), 0usize),
        (('C' as u8, 'G' as u8), 1usize),
        (('G' as u8, 'C' as u8), 1usize),
        (('C' as u8, 'T' as u8), 2usize),
        (('G' as u8, 'A' as u8), 2usize),
        (('T' as u8, 'A' as u8), 3usize),
        (('A' as u8, 'T' as u8), 3usize),
        (('T' as u8, 'C' as u8), 4usize),
        (('A' as u8, 'G' as u8), 4usize),
        (('T' as u8, 'G' as u8), 5usize),
        (('A' as u8, 'C' as u8), 5usize),
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
}

#[inline(always)]
fn is_insertion(a: Cigar) -> bool {
    return match a {
        Cigar::Ins(_) => true,
        _ => false,
    };
}
#[inline(always)]
fn query(a: Cigar) -> i64 {
    return match a {
        Cigar::Match(n) | Cigar::SoftClip(n) | Cigar::Ins(n) | Cigar::Diff(n) | Cigar::Equal(n) => {
            n as i64
        }
        _ => 0,
    };
}
#[inline(always)]
fn reference(a: Cigar) -> i64 {
    return match a {
        Cigar::Match(n) | Cigar::Del(n) | Cigar::Diff(n) | Cigar::Equal(n) | Cigar::RefSkip(n) => {
            n as i64
        }
        _ => 0,
    };
}

/// Return mapped parts of each read that overlap the other.
/// Returns A, B, genome coordiantes.
fn overlap_pieces(
    a: CigarStringView,
    b: CigarStringView,
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
                        },
                        Coordinates {
                            start: (b_read_pos + b_over) as u32,
                            stop: (b_read_pos + b_over + glen) as u32,
                        },
                        Coordinates {
                            start: genome_start as u32,
                            stop: genome_stop as u32,
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

    return result;
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{Cigar, CigarString};

    #[test]
    fn test_different_alignments() {
        let a = CigarString(vec![Cigar::Match(5), Cigar::Ins(3), Cigar::Match(5)]).into_view(0);
        let b = CigarString(vec![Cigar::Match(13)]).into_view(0);
        let r = overlap_pieces(a, b, false);
        dbg!(&r);
    }

    #[test]
    fn test_same_insertion() {
        let a = CigarString(vec![Cigar::Match(10), Cigar::Ins(8), Cigar::Match(10)]).into_view(8);
        let b = CigarString(vec![Cigar::Match(10), Cigar::Ins(8), Cigar::Match(10)]).into_view(5);
        let r = overlap_pieces(a, b, false);
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
        let r = overlap_pieces(a, b, true);
        let expected = [
            [
                Coordinates { start: 0, stop: 11 },
                Coordinates {
                    start: 10,
                    stop: 21,
                },
                Coordinates {
                    start: 10,
                    stop: 21,
                },
            ],
            [
                Coordinates {
                    start: 11,
                    stop: 23,
                },
                Coordinates {
                    start: 21,
                    stop: 33,
                },
                Coordinates {
                    start: 21,
                    stop: 33,
                },
            ],
            [
                Coordinates {
                    start: 23,
                    stop: 36,
                },
                Coordinates {
                    start: 33,
                    stop: 46,
                },
                Coordinates {
                    start: 33,
                    stop: 46,
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

        let r = overlap_pieces(a, b, true);

        let expected = [
            [
                Coordinates { start: 0, stop: 10 },
                Coordinates { start: 3, stop: 13 },
                Coordinates { start: 8, stop: 18 },
            ],
            [
                Coordinates {
                    start: 10,
                    stop: 67,
                },
                Coordinates {
                    start: 13,
                    stop: 70,
                },
                Coordinates {
                    start: 18,
                    stop: 75,
                },
            ],
            [
                Coordinates {
                    start: 67,
                    stop: 90,
                },
                Coordinates {
                    start: 70,
                    stop: 93,
                },
                Coordinates {
                    start: 75,
                    stop: 98,
                },
            ],
        ];

        assert_eq!(r, expected);
    }
}
