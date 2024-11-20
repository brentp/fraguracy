use rust_htslib::bam::{self, Read, Reader};
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::io::Result;
use std::path::PathBuf;

pub struct Denominator {
    heap: BinaryHeap<Reverse<u32>>,
    min_mapping_quality: u8,
    min_homopolymer_length: u8,
    max_homopolymer_distance: u8,
    contigs: Vec<String>,
    cur_tig: i32,
    cur_pos: u64,
    // tid, start, stop
    next_read: Option<(i32, usize, usize)>,
}

pub struct Depth {
    tid: i32,
    start: usize,
    stop: usize,
    depth: u32,
}

impl Denominator {
    pub fn new(
        contigs: Vec<String>,
        min_homopolymer_length: u8,
        max_homopolymer_distance: u8,
        min_mapping_quality: u8,
    ) -> Self {
        Self {
            heap: BinaryHeap::new(),
            min_mapping_quality,
            min_homopolymer_length,
            max_homopolymer_distance,
            contigs,
            cur_tig: 0,
            cur_pos: 0,
            next_read: None,
        }
    }

    pub fn add_read(&mut self, read: &bam::Record) -> Option<Depth> {
        if read.mapq() < self.min_mapping_quality {
            return None;
        }

        let tid = read.tid();
        let pos = read.pos();
        let mapq = read.mapq();
        None
    }
}

pub(crate) fn denominator_main(
    bam_path: PathBuf,
    min_homopolymer_length: u8,
    max_homopolymer_distance: u8,
    fasta: PathBuf,
    min_mapping_quality: u8,
    output_prefix: String,
) -> Result<()> {
    let mut bam = Reader::from_path(bam_path)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
    bam.set_reference(fasta)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
    bam.set_threads(1)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;

    let contigs = bam.header().target_names();

    Ok(())
}
