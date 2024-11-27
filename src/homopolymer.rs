pub(crate) const HP_REGEX: &str = "A{3,}|C{3,}|G{3,}|T{3,}";
use regex::Regex;
use rust_lapper::{Interval, Lapper};

/// Find the homopolymers and return an interval tree with the positions of the homopolymers.
pub(crate) fn find_homopolymers(seq: &[u8], re: &Regex) -> Lapper<u32, u8> {
    let seq_str = unsafe { std::str::from_utf8_unchecked(seq) };
    let matches = re.find_iter(seq_str);
    let intervals: Vec<Interval<u32, u8>> = matches
        .map(|m| Interval {
            start: m.range().start as u32,
            stop: m.range().end as u32,
            val: 0,
        })
        .collect();
    Lapper::new(intervals)
}

/// return a negative number if the hp is before the position, accounting for strand.
/// and 0 if the hp contains the position, otherwise a positive number.
///  hphphp---pos---->
///
pub(crate) fn hp_distance(
    hps: Option<&[&Interval<u32, u8>]>,
    pos: u32,
    read_start: u32,
    read_stop: u32,
    strand: i8,
) -> i8 {
    let mut dist = crate::fraguracy::MAX_HP_DIST;
    // strand will be 1 for forward, -1 for reverse
    for hp in hps.map(|hps| hps.iter()).unwrap_or_default() {
        // first we check if the hp is within 3 bases of the read start or stop.
        // since this could truncate the hp and not affect the read.
        // cases to exclude:
        // read: ----------->
        // hp:   AAAAA
        // pos:

        if hp.stop >= read_start && hp.stop < read_start + 3 {
            continue;
        }
        if hp.start < read_stop && hp.start > read_stop - 3 {
            continue;
        }

        assert!(pos >= read_start && pos <= read_stop);

        let d = if pos < hp.start {
            (hp.start - pos) as i64
        } else if pos > hp.stop {
            -((pos - hp.stop) as i64)
        } else {
            0i64
        };
        // now we check distance of pos to hp.
        let mut d = d.clamp(
            -crate::fraguracy::MAX_HP_DIST as i64,
            crate::fraguracy::MAX_HP_DIST as i64,
        ) as i8;
        if strand == -1 {
            d = -d;
        }
        //dbg!(pos, hp.start, hp.stop, strand, d);
        if d.abs() < dist.abs() {
            dist = d;
        }
    }
    /*
    if dist != crate::fraguracy::MAX_HP_DIST {
        eprintln!(
            "pos: {}, read: {}-{}, strand: {}, dist: {}",
            pos, read_start, read_stop, strand, dist
        );
    }
    */
    dist
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_homopolymers() {
        // Test sequence with various homopolymers
        let seq = b"AAATCCCGAAAGGGGTTTT";
        let re = Regex::new(HP_REGEX).expect("invalid regex");
        let homopolymers = find_homopolymers(seq, &re);

        // Convert results to vec for easier testing
        let results: Vec<_> = homopolymers.iter().collect();

        // Expected homopolymers: AAA, CCC, AAAA, GGGG, TTTT
        assert_eq!(results.len(), 5);

        // Check each homopolymer position
        assert_eq!(results[0].start, 0); // AAA
        assert_eq!(results[0].stop, 3);

        assert_eq!(results[1].start, 4); // CCC
        assert_eq!(results[1].stop, 7);

        assert_eq!(results[2].start, 8); // AAAA
        assert_eq!(results[2].stop, 11);

        assert_eq!(results[3].start, 11); // GGGG
        assert_eq!(results[3].stop, 15);

        assert_eq!(results[4].start, 15); // TTTT
        assert_eq!(results[4].stop, 19);
    }

    #[test]
    fn test_hp_distance() {
        let hp = vec![Interval {
            start: 9,
            stop: 12,
            val: 0,
        }];
        let hp_refs: Vec<&Interval<u32, u8>> = hp.iter().collect();

        // Test homopolymer near read end
        assert_eq!(
            hp_distance(Some(&hp_refs), 11, 10, 20, 1,),
            crate::fraguracy::MAX_HP_DIST
        );

        // Test forward strand
        assert_eq!(hp_distance(Some(&hp_refs), 15, 5, 20, 1,), -3);

        // Test reverse strand
        assert_eq!(hp_distance(Some(&hp_refs), 15, 5, 20, -1,), 3);
    }
}
