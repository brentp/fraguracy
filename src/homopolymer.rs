pub(crate) const HP_REGEX: &str = "A{3,}|C{3,}|G{3,}|T{3,}";
use coitrees::{Interval, IntervalTree};
use regex::Regex;

/// Find the homopolymers and return an interval tree with the positions of the homopolymers.
pub(crate) fn find_homopolymers(seq: &[u8], re: &Regex) -> coitrees::BasicCOITree<u32, u32> {
    let seq_str = unsafe { std::str::from_utf8_unchecked(seq) };
    let matches = re.find_iter(seq_str);
    let mut intervals: Vec<Interval<u32>> = matches
        .map(|m| Interval {
            first: m.range().start as i32,
            last: m.range().end as i32,
            metadata: 0,
        })
        .collect();

    // coitrees requires sorted intervals
    intervals.sort_by_key(|i| (i.first, i.last));

    log::info!("found {} homopolymers with regex: {re}", intervals.len());
    coitrees::BasicCOITree::new(&intervals)
}

/// return a negative number if the hp is before the position, accounting for strand.
/// and 0 if the hp contains the position, otherwise a positive number.
/// Returns None if the distance is greater than MAX_HP_DIST
///  hphphp---pos---->
///
pub(crate) fn hp_distance(
    hps: Option<&[coitrees::IntervalNode<u32, u32>]>,
    pos: u32,
    read_start: u32,
    read_stop: u32,
    _strand: i8,
) -> Option<i8> {
    let mut dist: Option<i8> = None;

    if let Some(hps) = hps {
        for hp in hps {
            // first we check if the hp is within 3 bases of the read start or stop
            /*
            if hp.last as u32 >= read_start && (hp.last as u32) < read_start + 3 {
                continue;
            }
            if hp.first as u32 <= read_stop && hp.first as u32 > read_stop - 3 {
                continue;
            }
            */

            assert!(pos >= read_start && pos <= read_stop);

            let d = if pos < hp.first as u32 {
                hp.first as i64 - pos as i64
            } else if pos > hp.last as u32 {
                -(pos as i64 - hp.last as i64)
            } else {
                0i64
            };

            if d < -crate::fraguracy::MAX_HP_DIST as i64 || d > crate::fraguracy::MAX_HP_DIST as i64
            {
                continue;
            }

            let d = d as i8;
            if dist.is_none() || d.abs() < dist.unwrap().abs() {
                dist = Some(d);
            }
        }
    }
    dist
}

#[cfg(test)]
mod tests {
    use super::*;
    use coitrees::IntervalNode;

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
        assert_eq!(results[0].first, 0); // AAA
        assert_eq!(results[0].last, 3);

        assert_eq!(results[1].first, 4); // CCC
        assert_eq!(results[1].last, 7);

        assert_eq!(results[2].first, 8); // AAAA
        assert_eq!(results[2].last, 11);

        assert_eq!(results[3].first, 11); // GGGG
        assert_eq!(results[3].last, 15);

        assert_eq!(results[4].first, 15); // TTTT
        assert_eq!(results[4].last, 19);
    }

    #[test]
    fn test_hp_distance() {
        let intervals = vec![IntervalNode::new(9, 12, 0)];
        let refs = intervals;

        // Test homopolymer near read end
        assert_eq!(hp_distance(Some(refs.as_slice()), 11, 10, 20, 1,), None);

        // Test forward strand
        assert_eq!(hp_distance(Some(refs.as_slice()), 15, 5, 20, 1,), Some(-3));

        // Test reverse strand
        assert_eq!(hp_distance(Some(refs.as_slice()), 15, 5, 20, -1,), Some(-3));

        // Test distant homopolymer
        assert_eq!(hp_distance(Some(refs.as_slice()), 115, 105, 120, 1,), None);
    }
}
