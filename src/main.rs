mod fraguracy;
#[macro_use]
extern crate lazy_static;
use clap::{Parser, Subcommand};
use std::path::PathBuf;

use rust_htslib::bam::{IndexedReader, Read, Reader};
use rust_htslib::faidx;
use rustc_hash::FxHashMap;

use crate::fraguracy::Stat;

use std::env;
use std::str;

#[derive(Debug, Parser)]
#[command(name = "fraguracy")]
#[command(about = "read accuracy from fragment overlaps")]

// see: https://docs.rs/clap/latest/clap/_cookbook/git_derive/index.html
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    #[command(arg_required_else_help = true)]
    Extract {
        bam: PathBuf,
        fasta: Option<PathBuf>,
        #[arg(short, long, default_value_t = 150)]
        max_read_length: u8,
        #[arg(short, long, default_value_t = 3)]
        bin_size: u8,
    },
}

fn main() {
    let args = Cli::parse();
    if env::var("RUST_LOG").is_err() {
        env::set_var("RUST_LOG", "info")
    }
    env_logger::init();

    match args.command {
        Commands::Extract {
            bam,
            fasta,
            bin_size,
            max_read_length,
        } => {
            extract_main(bam, fasta, bin_size as u32, max_read_length as u32);
        }
    }
}

fn extract_main(path: PathBuf, fasta_path: Option<PathBuf>, bin_size: u32, max_read_length: u32) {
    //let args: Vec<String> = env::args().collect();
    let mut map = FxHashMap::default();
    let min_base_qual = 10u8;
    let min_map_q = 10u8;

    let mut bam = Reader::from_path(&path).expect("error reading bam file {path}");
    bam.set_threads(3).expect("error setting threads");

    let mut ibam =
        IndexedReader::from_path(&path).expect("bam file (path) must be sorted and indexed");
    ibam.set_threads(3)
        .expect("error setting threads on indexed reader");

    let fasta: Option<faidx::Reader> = if let Some(fa_path) = fasta_path {
        bam.set_reference(&fa_path)
            .expect("Error setting reference for file");
        ibam.set_reference(&fa_path)
            .expect("Error setting reference for file");

        let fa = faidx::Reader::from_path(fa_path).expect("error opening faidx");
        Some(fa)
    } else {
        None
    };

    let bins = max_read_length / bin_size;
    let mut counts = fraguracy::Counts::new(ibam, bins as usize);

    let mut n_total = 0;
    let mut n_pairs = 0;
    let chroms: Vec<String> = bam
        .header()
        .target_names()
        .iter()
        .map(|n| unsafe { str::from_utf8_unchecked(n) }.to_string())
        .collect();

    let mut last_tid: i32 = 0;
    bam.rc_records()
        .map(|r| {
            n_total += 1;
            r.expect("error parsing read")
        })
        .filter(fraguracy::filter_read)
        .for_each(|b| {
            let name = unsafe { str::from_utf8_unchecked(b.qname()) }.to_string();
            if b.is_first_in_template() {
                n_pairs += 1;
            }
            if b.tid() != last_tid {
                log::info!("processed chromosome: {}", chroms[last_tid as usize]);
                last_tid = b.tid();
            }

            // by not checking the order here, we allow bams sorted by read name (with position flipped)
            // this gives about 5% performance penalty over checking b.pos() < b.mpos(), but allows us
            // to support more files.
            if !map.contains_key(&name) {
                assert!(!map.contains_key(&name));
                map.insert(name, b.clone());
            } else if let Some(a) = map.remove(&name) {
                if a.mapq() < min_map_q {
                    return;
                }
                if b.mapq() < min_map_q {
                    return;
                }
                // we know a is before b, but we don't know if they overlap.
                if a.cigar().end_pos() < b.pos() {
                    return;
                }
                let tid = a.tid() as usize;
                counts.increment(
                    a,
                    b,
                    min_base_qual,
                    min_map_q,
                    bin_size,
                    &fasta,
                    &chroms[tid],
                );
            } else {
                log::warn!("not found: {:?}{:?}", name, b.pos());
            }
        });
    log::info!(
        "[FINAL] map len:{:?} total reads: {:?}, pairs: {}, counts: {:?} \
        \n mismatches: {} matches: {}",
        map.len(),
        n_total,
        n_pairs,
        counts.cnts,
        counts.mismatches,
        counts.matches,
    );

    let stats = Stat::from_counts(counts, bin_size as usize);
    let header = Stat::header();
    println!("{header}");
    stats.iter().for_each(|s| println!("{s}"));
}
