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
}
