use rust_htslib::bam::{record::Cigar, IndexedReader, Read, Record};

struct CachedRead {
    pos: u32,
    end_pos: u32,
    /// True when CIGAR is a single M/=/X op, so `qpos = genome_pos - pos` with no CIGAR walk.
    simple_cigar: bool,
    rec: Record,
}

/// Forward-streaming pileup cache over an indexed BAM/CRAM.
/// Calling `ibam.fetch((tid, pos, pos+1))` on every mismatch repeatedly seeks and
/// re-decompresses the same 64 KiB BGZF blocks. Instead, we stream `ibam` linearly
/// forward across a small active window of reads and only seek on chromosome changes,
/// backward jumps, or forward jumps > 32 KiB.
pub(crate) struct PileupWindow {
    ibam: IndexedReader,
    cur_tid: i32,
    /// Lower bound on genomic positions for which all overlapping reads remain in `reads`.
    window_min_pos: u32,
    eof: bool,
    /// Next record already read from `ibam` (starts after the last queried `genome_pos`).
    peek: Option<Record>,
    /// Recycled `Record` buffer to avoid allocating a new htslib `bam1_t` for skipped reads.
    scratch: Record,
    reads: Vec<CachedRead>,
    /// Single-entry cache for repeated mismatch queries at the same hotspot coordinate.
    last_result: Option<(i32, u32, [u32; 5])>,
}

impl PileupWindow {
    pub(crate) fn new(ibam: IndexedReader) -> Self {
        Self {
            ibam,
            cur_tid: -1,
            window_min_pos: 0,
            eof: false,
            peek: None,
            scratch: Record::new(),
            reads: Vec::with_capacity(256),
            last_result: None,
        }
    }

    fn read_next(&mut self) {
        match self.ibam.read(&mut self.scratch) {
            Some(Ok(())) => {
                self.peek = Some(std::mem::replace(&mut self.scratch, Record::new()));
            }
            Some(Err(e)) => panic!("error reading bam in pileup window: {e}"),
            None => self.eof = true,
        }
    }

    /// Map 0-based reference `genome_pos` to 0-based query position in `r`, matching
    /// htslib's `Alignment::qpos()` (returns `None` inside deletions / ref skips).
    #[inline(always)]
    fn qpos_at_ref(r: &CachedRead, genome_pos: u32) -> Option<usize> {
        if r.simple_cigar {
            return Some((genome_pos - r.pos) as usize);
        }
        let mut rpos = r.pos;
        let mut qpos = 0usize;
        for op in r.rec.cigar().iter() {
            match *op {
                Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                    if genome_pos < rpos + n {
                        return Some(qpos + (genome_pos - rpos) as usize);
                    }
                    rpos += n;
                    qpos += n as usize;
                }
                Cigar::Ins(n) | Cigar::SoftClip(n) => qpos += n as usize,
                Cigar::Del(n) | Cigar::RefSkip(n) => {
                    if genome_pos < rpos + n {
                        return None;
                    }
                    rpos += n;
                }
                Cigar::HardClip(_) | Cigar::Pad(_) => {}
            }
        }
        None
    }

    /// Count high-quality `[A, C, G, T, N]` pileup bases at `(tid, genome_pos)`.
    pub(crate) fn query(
        &mut self,
        tid: i32,
        genome_pos: u32,
        min_map_qual: u8,
        min_base_qual: u8,
    ) -> [u32; 5] {
        if let Some((t, p, counts)) = self.last_result {
            if t == tid && p == genome_pos {
                return counts;
            }
        }

        // Reuse the open sequential stream if we haven't evicted reads covering `genome_pos`
        // and the next position is within 32 KiB (cheaper to scan forward than to BAI-seek).
        let can_stream = self.cur_tid == tid
            && genome_pos >= self.window_min_pos
            && (self.eof
                || self.peek.as_ref().is_some_and(|r| {
                    r.tid() == tid && (genome_pos as i64) <= r.pos() + 32_768
                }));

        if !can_stream {
            self.reads.clear();
            self.peek = None;
            self.eof = false;
            self.cur_tid = tid;
            self.window_min_pos = genome_pos;
            // Fetch to end of chromosome so subsequent forward queries stream without re-seeking.
            self.ibam
                .fetch((tid, genome_pos as i64, i64::MAX))
                .expect("Error seeking to genomic position");
            self.read_next();
        } else if genome_pos > self.window_min_pos + 1024 {
            // Advance `window_min_pos` and evict stale reads BEFORE pulling new ones so a
            // forward jump (up to 32 KiB) doesn't temporarily retain skipped intermediate reads.
            let new_min = genome_pos - 1024;
            self.window_min_pos = new_min;
            self.reads.retain(|r| r.end_pos > new_min);
        }

        // Pull in all reads starting at or before `genome_pos`.
        while let Some(rec) = self.peek.take() {
            if rec.tid() != tid || rec.pos() < 0 || (rec.pos() as u32) > genome_pos {
                self.peek = Some(rec);
                break;
            }
            // Match htslib default bam_plp mask: BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP.
            // Also require seq len >= 7 since bases within 3 bp of either read end are skipped below.
            if (rec.flags() & 0x704) == 0 && rec.mapq() >= min_map_qual && rec.qual().len() >= 7 {
                let cigar = rec.cigar();
                let end_pos = cigar.end_pos() as u32;
                if end_pos > self.window_min_pos {
                    let simple_cigar = cigar.len() == 1
                        && matches!(cigar[0], Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_));
                    let pos = rec.pos() as u32;
                    drop(cigar);
                    self.reads.push(CachedRead {
                        pos,
                        end_pos,
                        simple_cigar,
                        rec,
                    });
                } else {
                    self.scratch = rec;
                }
            } else {
                self.scratch = rec;
            }
            self.read_next();
        }

        let mut base_counts = [0u32; 5];
        for r in &self.reads {
            if r.pos > genome_pos || genome_pos >= r.end_pos {
                continue;
            }
            if let Some(qpos) = Self::qpos_at_ref(r, genome_pos) {
                let qual = r.rec.qual();
                // Skip bases within 3 bp of either read end.
                if qpos < 3 || qpos + 3 >= qual.len() || qual[qpos] < min_base_qual {
                    continue;
                }
                let base_idx = match unsafe { r.rec.seq().decoded_base_unchecked(qpos) } {
                    b'A' => 0,
                    b'C' => 1,
                    b'G' => 2,
                    b'T' => 3,
                    _ => 4,
                };
                base_counts[base_idx] += 1;
            }
        }
        self.last_result = Some((tid, genome_pos, base_counts));
        base_counts
    }
}
