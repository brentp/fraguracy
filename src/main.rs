use rust_htslib::bam::{
    record::{Cigar, CigarStringView},
    Read, Reader, Record,
};
use rustc_hash::FxHashMap;

use std::env;
use std::str;

fn main() {
    let args: Vec<String> = env::args().collect();
    let mut map = FxHashMap::default();

    let mut bam = Reader::from_path(&args[1]).expect("error reading bam file {args[1]}");
    bam.set_threads(3).expect("error setting threads");
    let mut n_total = 0;
    let mut n_pairs = 0;
    let mut n_overlapping = 0;
    let chroms: Vec<String> = bam
        .header()
        .target_names()
        .iter()
        .map(|n| unsafe { str::from_utf8_unchecked(n) }.to_string())
        .collect();

    bam.rc_records()
        .map(|r| r.expect("error parsing read"))
        .filter(|r| {
            n_total += 1;
            r.tid() == r.mtid()
                && r.tid() >= 0
                && !r.is_unmapped()
                && !r.is_mate_unmapped()
                && (r.pos() - r.mpos()).abs() < 1000
                && !r.is_supplementary()
                && !r.is_secondary()
        })
        .for_each(|b| {
            let name = unsafe { str::from_utf8_unchecked(b.qname()) }.to_string();
            if b.is_first_in_template() {
                n_pairs += 1;
            }

            if b.pos() < b.mpos() || (b.pos() == b.mpos() && !map.contains_key(&name)) {
                map.insert(name, b.clone());
            } else if let Some(a) = map.remove(&name) {
                // so we know a is before b, but we don't know if they overlap.
                if a.cigar().end_pos() < b.pos() {
                    return;
                }
                n_overlapping += 1;
                assert!(a.pos() <= b.pos());
                let pieces = overlap_pieces(a.cigar(), b.cigar(), false);
                if pieces.len() == 0 { return }
                let a_seq = a.seq();
                let b_seq = b.seq();
                let mut bases_overlap = 0;
                let mut mismatch_bases = 0;
                for [a_chunk, b_chunk] in pieces {
                    let mut bi = b_chunk.start as usize;
                    bases_overlap += b_chunk.stop - b_chunk.start;
                    for ai in a_chunk.start .. a_chunk.stop {
                        let a_base = unsafe { a_seq.decoded_base_unchecked(ai as usize) } ;
                        let b_base = unsafe { b_seq.decoded_base_unchecked(bi) } ;
                        mismatch_bases +=  (a_base != b_base) as i32;
                        bi += 1;
                    }

                }
                if mismatch_bases > 4 {
                    eprintln!(
                        "{qname} {a_tid}:{a_pos}-{a_end}({a_cigar}),Q:{a_qual} <-> {b_tid}:{b_pos}-{b_end}({b_cigar})Q:{b_qual}  mismatches:({mismatch_bases}/{bases_overlap})",
                        qname = str::from_utf8(a.qname()).unwrap(),
                        a_tid = chroms[a.tid() as usize],
                        a_pos = a.pos(),
                        a_end = a.cigar().end_pos(),
                        a_cigar = a.cigar(),
                        b_tid = chroms[b.tid() as usize],
                        b_pos = b.pos(),
                        b_end = b.cigar().end_pos(),
                        b_cigar = b.cigar(),
                        a_qual = a.mapq(),
                        b_qual = b.mapq(),
                    );
                    eprintln!(
                        "map len:{:?} total: {:?}, overlapping: {:?}",
                        map.len(),
                        n_pairs,
                        n_overlapping
                    );
                }
            }
        });
    eprintln!(
        "map len:{:?} total: {:?}, overlapping: {:?}, pairs: {}",
        map.len(),
        n_total,
        n_overlapping,
        n_pairs,
    );
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
