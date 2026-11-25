use crate::fraguracy::Q_LOOKUP;
use rust_htslib::{
    bgzf::{self, CompressionLevel},
    tpool::ThreadPool,
};
use std::collections::{BTreeMap, VecDeque};
use std::io::Write;

/// Per-genomic-position denominator depth bucket.
/// Base qualities are binned into 5 buckets (see `Counts::qual_to_bin`), so there are
/// 5 * 5 = 25 possible `(read_a_bq_bin, read_b_bq_bin)` pairs, stored flat at
/// `idx = aq * 5 + bq` to avoid allocating a `HashMap` per genomic coordinate.
/// Bit `idx` in `mask` is set iff `counts[idx] > 0` (25 bits fit in a `u32`),
/// allowing `flush_before` to skip empty slots in O(1) and visit only non-zero bins.
#[derive(Clone, Copy)]
struct DepthSlot {
    mask: u32,
    counts: [u32; 25],
}

impl DepthSlot {
    const ZERO: Self = Self {
        mask: 0,
        counts: [0; 25],
    };

    #[inline(always)]
    fn merge_from(&mut self, other: &DepthSlot) {
        let mut mask = other.mask;
        while mask != 0 {
            let idx = mask.trailing_zeros() as usize;
            mask &= mask - 1;
            self.counts[idx] += other.counts[idx];
        }
        self.mask |= other.mask;
    }
}

/// Cap the contiguous ring-buffer span so large coordinate jumps (e.g. spliced reads with
/// long introns, centromeres, or name-sorted BAMs) spill into the sparse `overflow` map
/// instead of allocating huge empty gaps or flushing prematurely out of order.
const MAX_DEPTH_WINDOW: u32 = 65_536;

struct LastDepthEntry {
    start_pos: u32,
    last_pos: u32,
    a_bin: u8,
    b_bin: u8,
    dp: u32,
}

/// Fast base-10 integer formatter into a byte buffer (avoids `core::fmt` overhead when
/// writing millions of denominator-depth BED lines). A `u32` has at most 10 digits.
#[inline(always)]
fn push_u32(buf: &mut Vec<u8>, mut v: u32) {
    let mut tmp = [0u8; 10];
    let mut i = 10;
    loop {
        i -= 1;
        tmp[i] = b'0' + (v % 10) as u8;
        v /= 10;
        if v == 0 {
            break;
        }
    }
    buf.extend_from_slice(&tmp[i..]);
}

/// Sliding-window denominator-depth accumulator and run-length-encoded BGZF BED writer.
/// `slots[i]` corresponds to 0-based genomic coordinate `start_pos + i`.
/// Coordinates outside `MAX_DEPTH_WINDOW` of `start_pos` (such as downstream exons of spliced
/// reads spanning >65 kb introns) spill into `overflow` and are merge-drained in sorted order
/// by `flush_before`, keeping output strictly coordinate-sorted with zero overhead when unused.
/// Buffering directly into `out_buf` (64 KiB) avoids a redundant `BufWriter` memcpy per line.
pub(crate) struct DepthTracker {
    start_pos: u32,
    slots: VecDeque<DepthSlot>,
    overflow: BTreeMap<u32, DepthSlot>,
    last_entry: Option<LastDepthEntry>,
    writer: bgzf::Writer,
    out_buf: Vec<u8>,
}

impl DepthTracker {
    pub(crate) fn new(path: &str) -> std::io::Result<Self> {
        let mut writer = bgzf::Writer::from_path_with_level(path, CompressionLevel::Level(1))
            .map_err(std::io::Error::other)?;
        // `bgzf::Writer::set_thread_pool` clones and retains its own `ThreadPool` handle.
        if let Ok(tpool) = ThreadPool::new(2) {
            let _ = writer.set_thread_pool(&tpool);
        }
        let mut out_buf = Vec::with_capacity(64 * 1024 + 128);
        out_buf
            .extend_from_slice(b"#chrom\tstart\tend\tread1_bq_bin\tread2_bq_bin\tpair-ovl-depth\n");
        Ok(Self {
            start_pos: 0,
            slots: VecDeque::with_capacity(512),
            overflow: BTreeMap::new(),
            last_entry: None,
            writer,
            out_buf,
        })
    }

    /// Return a mutable slot for `pos` if it falls within `MAX_DEPTH_WINDOW` of `start_pos`,
    /// extending the front or back of `slots` as needed; otherwise `None` (caller uses `overflow`).
    #[inline(always)]
    fn slot_mut(&mut self, pos: u32) -> Option<&mut DepthSlot> {
        if self.slots.is_empty() {
            self.start_pos = pos;
            self.slots.push_back(DepthSlot::ZERO);
            return Some(&mut self.slots[0]);
        }
        if pos >= self.start_pos {
            let off = pos - self.start_pos;
            if off < MAX_DEPTH_WINDOW {
                let off = off as usize;
                if off >= self.slots.len() {
                    self.slots.resize(off + 1, DepthSlot::ZERO);
                }
                return Some(&mut self.slots[off]);
            }
        } else {
            let delta = self.start_pos - pos;
            if delta < MAX_DEPTH_WINDOW {
                for _ in 0..delta {
                    self.slots.push_front(DepthSlot::ZERO);
                }
                self.start_pos = pos;
                return Some(&mut self.slots[0]);
            }
        }
        None
    }

    #[inline(always)]
    fn append_entry(&mut self, bchrom: &str, entry: &LastDepthEntry) {
        self.out_buf.extend_from_slice(bchrom.as_bytes());
        self.out_buf.push(b'\t');
        push_u32(&mut self.out_buf, entry.start_pos);
        self.out_buf.push(b'\t');
        push_u32(&mut self.out_buf, entry.last_pos + 1);
        self.out_buf.push(b'\t');
        self.out_buf
            .extend_from_slice(Q_LOOKUP[entry.a_bin as usize].as_bytes());
        self.out_buf.push(b'\t');
        self.out_buf
            .extend_from_slice(Q_LOOKUP[entry.b_bin as usize].as_bytes());
        self.out_buf.push(b'\t');
        push_u32(&mut self.out_buf, entry.dp);
        self.out_buf.push(b'\n');
        if self.out_buf.len() >= 64 * 1024 {
            self.writer
                .write_all(&self.out_buf)
                .expect("error writing to bgzf file");
            self.out_buf.clear();
        }
    }

    /// Record one overlapping read pair at genomic coordinate `pos` with quantized
    /// quality bins `aq` and `bq` (each in `0..5`).
    #[inline(always)]
    pub(crate) fn add(&mut self, pos: u32, aq: u8, bq: u8) {
        // Flatten the 5x5 `(aq, bq)` matrix into `0..25` (row-major; decoded via `idx / 5, idx % 5`).
        let idx = aq as usize * 5 + bq as usize;
        // Positions outside `MAX_DEPTH_WINDOW` spill into `overflow` instead of flushing early.
        let slot = match self.slot_mut(pos) {
            Some(s) => s,
            None => self.overflow.entry(pos).or_insert(DepthSlot::ZERO),
        };
        slot.counts[idx] += 1;
        // Mark bit `idx` so `flush_before` visits this `(aq, bq)` bin without scanning all 25.
        slot.mask |= 1u32 << idx;
    }

    #[inline(always)]
    fn emit_slot(&mut self, bchrom: &str, pos: u32, slot: DepthSlot) {
        let mut mask = slot.mask;
        // Visit only set bits of `mask` in ascending `idx` order:
        // `trailing_zeros()` gives the lowest set bit; `mask &= mask - 1` clears it.
        while mask != 0 {
            let idx = mask.trailing_zeros() as usize;
            mask &= mask - 1;
            let dp = slot.counts[idx];
            let (aq, bq) = ((idx / 5) as u8, (idx % 5) as u8);

            if let Some(ref mut last) = self.last_entry {
                if aq == last.a_bin && bq == last.b_bin && dp == last.dp && last.last_pos + 1 == pos
                {
                    last.last_pos = pos;
                    continue;
                }
                let prev = std::mem::replace(
                    last,
                    LastDepthEntry {
                        start_pos: pos,
                        last_pos: pos,
                        a_bin: aq,
                        b_bin: bq,
                        dp,
                    },
                );
                self.append_entry(bchrom, &prev);
            } else {
                self.last_entry = Some(LastDepthEntry {
                    start_pos: pos,
                    last_pos: pos,
                    a_bin: aq,
                    b_bin: bq,
                    dp,
                });
            }
        }
    }

    /// Pop and run-length-merge all depth slots strictly before `bpos` (`i64::MAX` flushes all).
    pub(crate) fn flush_before(&mut self, bchrom: &str, bpos: i64) {
        loop {
            let ring =
                (!self.slots.is_empty() && (self.start_pos as i64) < bpos).then_some(self.start_pos);
            let ov = self
                .overflow
                .first_key_value()
                .map(|(&p, _)| p)
                .filter(|&p| (p as i64) < bpos);
            let Some(pos) = ring.into_iter().chain(ov).min() else {
                break;
            };
            let mut slot = DepthSlot::ZERO;
            if ring == Some(pos) {
                slot = self.slots.pop_front().unwrap();
                self.start_pos = pos.wrapping_add(1);
            }
            if ov == Some(pos) {
                slot.merge_from(&self.overflow.remove(&pos).unwrap());
            }
            if slot.mask != 0 {
                self.emit_slot(bchrom, pos, slot);
            }
        }
        if bpos == i64::MAX {
            if let Some(last) = self.last_entry.take() {
                self.append_entry(bchrom, &last);
            }
            if !self.out_buf.is_empty() {
                self.writer
                    .write_all(&self.out_buf)
                    .expect("error writing to bgzf file");
                self.out_buf.clear();
            }
        }
    }
}

impl Drop for DepthTracker {
    fn drop(&mut self) {
        if !self.out_buf.is_empty() {
            let _ = self.writer.write_all(&self.out_buf);
            self.out_buf.clear();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    #[test]
    fn test_spliced_reads_stay_sorted_and_merged() {
        let tmp = std::env::temp_dir().join(format!(
            "fraguracy-spliced-depth-test-{}.bed.gz",
            std::process::id()
        ));
        let tmp_str = tmp.to_string_lossy().into_owned();

        let mut dt = DepthTracker::new(&tmp_str).unwrap();

        // Pair 1 (spliced across a 100 kb intron > MAX_DEPTH_WINDOW): exons [1000, 1050) and [101050, 101100)
        for p in 1000..1050 {
            dt.add(p, 2, 2);
        }
        for p in 101050..101100 {
            dt.add(p, 2, 2);
        }
        dt.flush_before("chr1", 1000);

        // Pair 2 (normal read overlapping exon 1): [1010, 1060)
        for p in 1010..1060 {
            dt.add(p, 2, 2);
        }
        dt.flush_before("chr1", 1010);

        // Pair 3 (normal read overlapping exon 2 after a 100 kb gap): [101060, 101110)
        for p in 101060..101110 {
            dt.add(p, 2, 2);
        }
        dt.flush_before("chr1", 101060);

        // Pair 4 (after a 500 kb centromere jump): [600000, 600010)
        for p in 600000..600010 {
            dt.add(p, 3, 3);
        }
        dt.flush_before("chr1", 600000);
        dt.flush_before("chr1", i64::MAX);
        drop(dt);

        let mut gz = flate2::read::MultiGzDecoder::new(std::fs::File::open(&tmp).unwrap());
        let mut bed = String::new();
        gz.read_to_string(&mut bed).unwrap();
        let _ = std::fs::remove_file(&tmp);

        let rows: Vec<&str> = bed.lines().collect();
        assert_eq!(
            rows,
            vec![
                "#chrom\tstart\tend\tread1_bq_bin\tread2_bq_bin\tpair-ovl-depth",
                "chr1\t1000\t1010\t20-36\t20-36\t1",
                "chr1\t1010\t1050\t20-36\t20-36\t2",
                "chr1\t1050\t1060\t20-36\t20-36\t1",
                "chr1\t101050\t101060\t20-36\t20-36\t1",
                "chr1\t101060\t101100\t20-36\t20-36\t2",
                "chr1\t101100\t101110\t20-36\t20-36\t1",
                "chr1\t600000\t600010\t37-59\t37-59\t1",
            ]
        );
    }
}
