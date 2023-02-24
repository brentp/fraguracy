use ndarray::prelude::{Array, ArrayBase, ArrayView6};
use ndarray::{Array6, ArrayViewMut6};
use rust_htslib::bam::{
    record::{Cigar, CigarStringView},
    Read, Reader, Record,
};
use std::collections::HashMap;
use std::rc::Rc;
use std::str;

pub(crate) struct Counts {
    //  read, pos, mq, bp, ctx{6} */
    pub muts: Array6<u64>,
    //  read, pos, mq, bp, ctx{2} */
    pub cnts: Array6<u64>,
}

impl Counts {
    pub(crate) fn new() -> Self {
        Counts {
            /*                         read1/2, F/R, pos, mq, bq, ctx */
            cnts: Array::zeros((2, 2, 50, 5, 5, 2)),
            muts: Array::zeros((2, 2, 50, 5, 5, 6)),
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

    pub(crate) fn increment(&mut self, a: Rc<Record>, b: Rc<Record>, min_base_qual: u8) {
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

        for [a_chunk, b_chunk] in pieces {
            let mut bi = b_chunk.start as usize;
            for ai in a_chunk.start..a_chunk.stop {
                let aq = a_qual[ai as usize];
                if aq < min_base_qual {
                    bi += 1;
                    continue;
                }
                let bq = b_qual[bi as usize];
                if bq < min_base_qual {
                    bi += 1;
                    continue;
                }
                let bq = Counts::qual_to_bin(bq);
                let aq = Counts::qual_to_bin(aq);

                let a_base = unsafe { a_seq.decoded_base_unchecked(ai as usize) };
                let b_base = unsafe { b_seq.decoded_base_unchecked(bi) };

                let mismatch = a_base != b_base;

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

                if a_base != b_base {
                    // TODO: pileup and vote to determine error.
                    let index = [
                        a.is_first_in_template() as usize,
                        1 - (a.is_reverse() as usize),
                        a_pos,
                        amq as usize,
                        aq as usize,
                    ];
                }
            }
        }
        /*
        if mismatch_bases > 0 {
            eprintln!(
                "{qname} {a_tid}:{a_pos}-{a_end}({a_cigar}),Q:{a_qual} \
                        <-> {b_tid}:{b_pos}-{b_end}({b_cigar})Q:{b_qual} \
                          mismatches:({mismatch_bases}/{bases_overlap}) quals: {mquals:?}",
                qname = str::from_utf8(a.qname()).unwrap(),
                mismatch_bases = mismatch_bases as u64,
                a_tid = a.tid(),
                a_pos = a.pos(),
                a_end = a.cigar().end_pos(),
                a_cigar = a.cigar(),
                b_tid = b.tid(),
                b_pos = b.pos(),
                b_end = b.cigar().end_pos(),
                b_cigar = b.cigar(),
                a_qual = a.mapq(),
                b_qual = b.mapq(),
                mquals = mquals,
            );
            eprintln!(
                "map len:{:?} total: {:?}, overlapping-bases: {:?}",
                map.len(),
                n_pairs,
                bases_overlapping,
            );
        }
            */
    }
}

lazy_static! {
    pub(crate) static ref CONTEXT_LOOKUP: HashMap<(u8, u8), u8> = HashMap::from([
        (('C' as u8, 'A' as u8), 0u8),
        (('G' as u8, 'T' as u8), 0u8),
        (('C' as u8, 'G' as u8), 1u8),
        (('G' as u8, 'C' as u8), 1u8),
        (('C' as u8, 'T' as u8), 2u8),
        (('G' as u8, 'A' as u8), 2u8),
        (('T' as u8, 'A' as u8), 3u8),
        (('A' as u8, 'T' as u8), 3u8),
        (('T' as u8, 'C' as u8), 4u8),
        (('A' as u8, 'G' as u8), 4u8),
        (('T' as u8, 'G' as u8), 5u8),
        (('A' as u8, 'C' as u8), 5u8),
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

fn b() {
    let mut m: Array6<u64> = Array::zeros((2, 2, 50, 5, 5, 6));

    start(m.view_mut());
    start(m.view_mut());
}

fn start(mut counts: ArrayViewMut6<u64>) {
    //let mut counts = Array::zeros((2, 2, 2, 50, 5, 5));
    counts[[0, 0, 0, 20, 2, 2]] += 1;
}

struct ReadInfo {
    read_pos: u8,
    base_q: u8,
    map_q: u8,
    context: u8,
    read: u8, // read 1 or read 2.
    count: u64,
}

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct ReadCoordinates {
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

/// Return mapped parts of each read that overlap the other. Coordinates are in read-space.
fn overlap_pieces(
    a: CigarStringView,
    b: CigarStringView,
    skip_insertions: bool,
) -> Vec<[ReadCoordinates; 2]> {
    let aend = a.end_pos();
    let bend = b.end_pos();
    //let astart = a.pos() + a.leading_softclips();
    //let bstart = b.pos() + b.leading_softclips();
    if aend <= b.pos() || bend <= a.pos() {
        return vec![];
    }
    let mut result: Vec<[ReadCoordinates; 2]> = Vec::new();
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
                        ReadCoordinates {
                            start: (a_read_pos + a_over) as u32,
                            stop: (a_read_pos + a_over + glen) as u32,
                        },
                        ReadCoordinates {
                            start: (b_read_pos + b_over) as u32,
                            stop: (b_read_pos + b_over + glen) as u32,
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
                ReadCoordinates { start: 0, stop: 11 },
                ReadCoordinates {
                    start: 10,
                    stop: 21,
                },
            ],
            [
                ReadCoordinates {
                    start: 11,
                    stop: 23,
                },
                ReadCoordinates {
                    start: 21,
                    stop: 33,
                },
            ],
            [
                ReadCoordinates {
                    start: 23,
                    stop: 36,
                },
                ReadCoordinates {
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
                ReadCoordinates { start: 0, stop: 10 },
                ReadCoordinates { start: 3, stop: 13 },
            ],
            [
                ReadCoordinates {
                    start: 10,
                    stop: 67,
                },
                ReadCoordinates {
                    start: 13,
                    stop: 70,
                },
            ],
            [
                ReadCoordinates {
                    start: 67,
                    stop: 90,
                },
                ReadCoordinates {
                    start: 70,
                    stop: 93,
                },
            ],
        ];

        assert_eq!(r, expected);
    }
}
