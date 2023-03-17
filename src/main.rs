mod files;
mod fraguracy;
//mod plot;
#[macro_use]
extern crate lazy_static;
use clap::{Parser, Subcommand};
use linear_map::LinearMap;

use std::collections::HashMap;
use std::path::PathBuf;

use rust_htslib::bam;
use rust_htslib::bam::{IndexedReader, Read, Reader};
use rust_htslib::faidx;
use rustc_hash::FxHashMap;

use rayon::prelude::*;

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
        #[arg(
            short,
            long,
            help = "fasta for use with crams and/or to use as 'truth'"
        )]
        fasta: Option<PathBuf>,
        bams: Vec<PathBuf>,
        #[arg(
            short,
            long,
            default_value_t = String::from("fraguracy-"),
            help = "prefix for output files"
        )]
        output_prefix: String,

        #[arg(
            short,
            long,
            default_value_t = 150,
            help = "indicate the maximum read length in the alignment file"
        )]
        max_read_length: u8,
        #[arg(
            short,
            long,
            default_value_t = 3,
            help = "parition the read into chunks/bins of this size"
        )]
        bin_size: u8,
        #[arg(
            short = 'Q',
            long,
            default_value_t = 50,
            help = "only consider pairs where both reads have this mapping-quality or higher (good to leave this high)"
        )]
        min_mapping_quality: u8,
    },
    //Plot { tsv: PathBuf, },
}

fn get_sample_name(hmap: HashMap<String, Vec<LinearMap<String, String>>>) -> String {
    if let Some(lm) = hmap.get("RG") {
        let sm = String::from("SM");
        if let Some(v) = lm[0].get(&sm) {
            return (*v).clone();
        } else {
            return String::from("");
        }
    } else {
        return String::from("");
    }
}

fn main() {
    let args = Cli::parse();
    if env::var("RUST_LOG").is_err() {
        env::set_var("RUST_LOG", "info")
    }
    env_logger::init();

    match args.command {
        Commands::Extract {
            bams,
            fasta,
            output_prefix,
            bin_size,
            max_read_length,
            min_mapping_quality,
        } => {
            extract_main(
                bams,
                fasta,
                PathBuf::from(output_prefix),
                bin_size as u32,
                max_read_length as u32,
                min_mapping_quality,
            );
        } //Commands::Plot { tsv } => plot::plot(tsv),
    }
}

fn process_bam(
    path: PathBuf,
    fasta_path: Option<PathBuf>,
    output_prefix: PathBuf,
    bin_size: u32,
    max_read_length: u32,
    min_mapping_quality: u8,
    min_base_qual: u8,
) {
    let mut bam = Reader::from_path(&path).expect("error reading bam file {path}");
    bam.set_threads(3).expect("error setting threads");
    let mut map = FxHashMap::default();

    let mut ibam =
        IndexedReader::from_path(&path).expect("bam file (path) must be sorted and indexed");
    ibam.set_threads(3)
        .expect("error setting threads on indexed reader");

    let fasta: Option<faidx::Reader> = if let Some(fa_path) = fasta_path.clone() {
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

    let hmap = bam::Header::from_template(bam.header()).to_hashmap();
    let sample_name = get_sample_name(hmap);
    log::info!("found sample {sample_name}");
    let output_prefix: PathBuf =
        (output_prefix.to_string_lossy().to_string() + &sample_name).into();

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
                if a.mapq() < min_mapping_quality {
                    return;
                }
                if b.mapq() < min_mapping_quality {
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
                    min_mapping_quality,
                    bin_size,
                    &fasta,
                    &chroms[tid],
                );
            } else {
                log::warn!("not found: {:?}{:?}", name, b.pos());
            }
        });
    log::info!(
        "[FINAL] map len:{:?} total reads: {:?}, pairs: {} \
        \n mismatches: {} matches: {}",
        map.len(),
        n_total,
        n_pairs,
        counts.mismatches,
        counts.matches,
    );

    let stats = Stat::from_counts(&counts, bin_size as usize);
    files::write_stats(stats, output_prefix.clone());
    files::write_errors(&counts, output_prefix, chroms);
}

fn extract_main(
    paths: Vec<PathBuf>,
    fasta_path: Option<PathBuf>,
    output_prefix: PathBuf,
    bin_size: u32,
    max_read_length: u32,
    min_mapping_quality: u8,
) {
    //let args: Vec<String> = env::args().collect();
    let min_base_qual = 10u8;

    paths.par_iter().for_each(|path| {
        process_bam(
            path.clone(),
            fasta_path.clone(),
            output_prefix.clone(),
            bin_size,
            max_read_length,
            min_mapping_quality,
            min_base_qual,
        )
    })
}
