mod combine_counts;
mod combine_errors;
mod homopolymer;

mod files;
mod fraguracy;

//mod plot;
#[macro_use]
extern crate lazy_static;
use clap::{Parser, Subcommand};
use fraguracy::ConfidenceInterval;
use homopolymer::find_homopolymers;
use linear_map::LinearMap;
use regex::Regex;

use rust_lapper::{Interval, Lapper};

use std::io::BufRead;

type Iv = Interval<u32, u32>;

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
#[command(
    version,
    about = "fraguracy: unbiased error profile analysis for short read sequencing",
    author = "Brent S Pedersen",
    help_template = "{about}\nversion:{version}\n\n{usage-heading} {usage} \n\nOPTIONS:\n{options}\n\n\x1b[1m\x1b[4mCOMMANDS:\x1b[0m\n{subcommands}"
)]
#[command(arg_required_else_help = true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    #[command(
        arg_required_else_help = true,
        about = "combine error bed files from extract"
    )]
    CombineErrors {
        #[arg(
            short,
            long,
            required = true,
            help = "path for to fai (not fasta) file"
        )]
        fai_path: PathBuf,

        #[arg(help = "path to error bed files from extract")]
        errors: Vec<PathBuf>,

        #[arg(
            short,
            long,
            default_value_t = String::from("fraguracy-combined-errors.bed"),
            help = "path for output bed file"
        )]
        output_path: String,
    },

    #[command(
        arg_required_else_help = true,
        about = "combine counts.txt files from extract"
    )]
    CombineCounts {
        #[arg(help = "path to counts.txt files from extract")]
        counts: Vec<PathBuf>,

        #[arg(
            short,
            long,
            default_value_t = String::from("fraguracy-combined-counts.txt"),
            help = "path for output counts file"
        )]
        output_path: String,
    },

    #[command(
        arg_required_else_help = true,
        about = "error profile pair overlaps in bam/cram"
    )]
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
            help = "restrict analysis to the regions given in this BED file"
        )]
        regions: Option<PathBuf>,

        #[arg(
            short,
            long,
            help = "exclude from analysis the regions given in this BED file"
        )]
        exclude_regions: Option<PathBuf>,

        #[arg(
            short,
            long,
            default_value_t = 151,
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

        #[arg(
            short,
            long = "ci",
            help = "method for confidence interval calculation (see rust bpci crate)",
            default_value = "agresti-coull"
        )]
        ci: ConfidenceInterval,

        #[arg(
            short,
            long,
            help = "do not calculate denominator. This can shorten runtime.",
            default_value_t = false
        )]
        no_denominator: bool,

        #[arg(
            short = 'H',
            long,
            help = format!(
                "regex for homopolymer sequence to consider if denominator is calculated[default: {}]",
                homopolymer::HP_REGEX
            ),
            default_value = homopolymer::HP_REGEX
        )]
        homopolymer_regex: String,

        #[arg(
            short = 't',
            long,
            help = "use reference base as 'truth'",
            default_value_t = false
        )]
        reference_as_truth: bool,
    },
    //Plot { tsv: PathBuf, },
}

fn get_sample_name(hmap: HashMap<String, Vec<LinearMap<String, String>>>) -> String {
    if let Some(lm) = hmap.get("RG") {
        let sm = String::from("SM");
        if let Some(v) = lm[0].get(&sm) {
            (*v).clone()
        } else {
            String::from("")
        }
    } else {
        String::from("")
    }
}

fn main() -> std::io::Result<()> {
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
            regions,
            exclude_regions,
            bin_size,
            max_read_length,
            min_mapping_quality,
            ci,
            reference_as_truth,
            no_denominator,
            homopolymer_regex,
        } => extract_main(
            bams,
            fasta,
            PathBuf::from(output_prefix),
            regions,
            exclude_regions,
            bin_size as u32,
            max_read_length as u32,
            min_mapping_quality,
            ci,
            reference_as_truth,
            no_denominator,
            homopolymer_regex,
        ), //Commands::Plot { tsv } => plot::plot(tsv),
        Commands::CombineErrors {
            fai_path,
            errors,
            output_path,
        } => combine_errors::combine_errors_main(errors, fai_path, output_path),

        Commands::CombineCounts {
            counts,
            output_path,
        } => combine_counts::combine_counts_main(counts, output_path),
    }
}

fn read_bed(path: Option<PathBuf>) -> Option<HashMap<String, Lapper<u32, u32>>> {
    path.as_ref()?;

    let reader = files::open_file(path);
    reader.as_ref()?;
    let mut bed = HashMap::new();

    reader
        .expect("checked that reader is available")
        .lines()
        .for_each(|l| {
            let line = l.expect("error reading line");
            let fields: Vec<_> = line.split('\t').collect();
            if let (Ok(start), Ok(stop)) = (fields[1].parse::<u32>(), fields[2].parse::<u32>()) {
                let iv = Iv {
                    start,
                    stop,
                    val: 0,
                };
                let chrom = String::from(fields[0]);
                bed.entry(chrom).or_insert(Vec::new()).push(iv);
            }
        });

    let mut tree: HashMap<String, Lapper<u32, u32>> = HashMap::new();

    for (chrom, ivs) in bed.iter() {
        let ivs = ivs.clone();
        let chrom = (*chrom).clone();
        tree.insert(chrom, Lapper::new(ivs));
    }
    Some(tree)
}

fn get_tree<'a>(
    regions: &'a Option<HashMap<String, Lapper<u32, u32>>>,
    chrom: &String,
) -> Option<&'a Lapper<u32, u32>> {
    let tree: Option<&Lapper<u32, u32>> = if let Some(r) = regions {
        r.get(chrom)
    } else {
        None
    };
    tree
}

#[allow(clippy::too_many_arguments)]
fn process_bam(
    path: PathBuf,
    fasta_path: Option<PathBuf>,
    regions: Option<PathBuf>,
    exclude_regions: Option<PathBuf>,
    bin_size: u32,
    max_read_length: u32,
    min_mapping_quality: u8,
    min_base_qual: u8,
    reference_as_truth: bool,
    output_prefix: PathBuf,
    no_denominator: bool,
    homopolymer_regex: Option<Regex>,
) -> (fraguracy::InnerCounts, Vec<String>, String) {
    let mut bam =
        Reader::from_path(&path).unwrap_or_else(|_| panic!("error reading bam file {path:?}"));
    bam.set_threads(1).expect("error setting threads");
    let mut map = FxHashMap::default();

    let include_regions = read_bed(regions);
    let exclude_regions = read_bed(exclude_regions);

    let mut ibam = IndexedReader::from_path(&path)
        .unwrap_or_else(|_| panic!("bam file {path:?} must be sorted and indexed"));
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

    let bins = (max_read_length as f64 / bin_size as f64).ceil() as u32;
    let mut counts = fraguracy::Counts::new(
        if reference_as_truth { None } else { Some(ibam) },
        bins as usize,
    );
    let hmap = bam::Header::from_template(bam.header()).to_hashmap();
    let sample_name = get_sample_name(hmap);
    log::info!("found sample {sample_name}");

    if !no_denominator {
        counts
            .set_depth_writer(
                &(output_prefix.to_string_lossy().to_string()
                    + &sample_name
                    + "-fraguracy-denominator-depth.bed.gz")
                    .to_string(),
            )
            .unwrap();
    }

    let mut n_total = 0;
    let mut n_pairs = 0;
    let chroms: Vec<String> = bam
        .header()
        .target_names()
        .iter()
        .map(|n| unsafe { str::from_utf8_unchecked(n) }.to_string())
        .collect();

    let mut include_tree: Option<&Lapper<u32, u32>> = get_tree(&include_regions, &chroms[0]);
    let mut exclude_tree: Option<&Lapper<u32, u32>> = get_tree(&exclude_regions, &chroms[0]);
    let mut hp_tree: Option<Lapper<u32, u8>> = if let Some(ref re) = homopolymer_regex {
        let chrom_seq = fasta
            .as_ref()
            .unwrap()
            .fetch_seq(&chroms[0], 0, i64::MAX as usize)
            .expect("error fetching sequence");
        Some(find_homopolymers(&chrom_seq, re))
    } else {
        None
    };

    let mut last_tid: i32 = -1;
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
                if last_tid != -1 {
                    log::info!(
                        "processed chromosome: {} unprocessed orphan pairs: {}",
                        chroms[last_tid as usize],
                        map.len()
                    );
                }
                last_tid = b.tid();

                // process the remaining entries in the hashmap in last_depth.
                counts.handle_depth(&chroms[last_tid as usize], i64::MAX);

                if include_regions.is_some() {
                    include_tree = get_tree(&include_regions, &chroms[last_tid as usize]);
                }
                if exclude_regions.is_some() {
                    exclude_tree = get_tree(&exclude_regions, &chroms[last_tid as usize]);
                }

                if hp_tree.is_some() {
                    let chrom_seq = fasta
                        .as_ref()
                        .unwrap()
                        .fetch_seq(&chroms[last_tid as usize], 0, i64::MAX as usize)
                        .expect("error fetching sequence from fasta.");
                    hp_tree = Some(find_homopolymers(
                        &chrom_seq,
                        homopolymer_regex.as_ref().unwrap(),
                    ));
                }
            }

            // by not checking the order here, we allow bams sorted by read name (with position flipped)
            // this gives about 5% performance penalty over checking b.pos() < b.mpos(), but allows us
            // to support more files.
            match map.entry(name) {
                std::collections::hash_map::Entry::Vacant(e) => {
                    e.insert(b);
                }
                std::collections::hash_map::Entry::Occupied(e) => {
                    let a = e.remove();

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
                        &include_tree,
                        &exclude_tree,
                        &hp_tree,
                    );
                }
            }
        });
    log::info!(
        "[FINAL] map len:{:?} total reads: {:?}, pairs: {} \
        \n mismatches: {} matches: {}",
        map.len(),
        n_total,
        n_pairs,
        counts.counts.mismatches,
        counts.counts.matches,
    );

    (counts.counts, chroms, sample_name)
}

#[allow(clippy::too_many_arguments)]
fn extract_main(
    paths: Vec<PathBuf>,
    fasta_path: Option<PathBuf>,
    output_prefix: PathBuf,
    regions: Option<PathBuf>,
    exclude_regions: Option<PathBuf>,
    bin_size: u32,
    max_read_length: u32,
    min_mapping_quality: u8,
    ci: ConfidenceInterval,
    reference_as_truth: bool,
    no_denominator: bool,
    homopolymer_regex: String,
) -> std::io::Result<()> {
    //let args: Vec<String> = env::args().collect();
    let min_base_qual = 20u8;

    let mut homopolymer_regex =
        Some(Regex::new(&homopolymer_regex).expect("error compiling homopolymer regex"));

    if no_denominator {
        homopolymer_regex = None;
    } else if fasta_path.is_none() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "fasta path must be provided if denominator is calculated",
        ));
    }

    let total_counts = paths
        .par_iter()
        .map(|path| {
            let (c, chroms, sample_name) = process_bam(
                path.clone(),
                fasta_path.clone(),
                regions.clone(),
                exclude_regions.clone(),
                bin_size,
                max_read_length,
                min_mapping_quality,
                min_base_qual,
                reference_as_truth,
                output_prefix.clone(),
                no_denominator,
                homopolymer_regex.clone(),
            );
            let output_prefix: PathBuf =
                (output_prefix.to_string_lossy().to_string() + &sample_name + "-").into();

            let stats = Stat::from_counts(&c, bin_size as usize, ci.clone());
            files::write_stats(stats, output_prefix.clone());
            files::write_errors(&c, output_prefix, chroms);

            c
        })
        .reduce_with(|mut a, b| {
            a += b;
            a
        });

    let total_counts = total_counts.expect("error accumulating total counts");
    if paths.len() > 1 {
        let bam = Reader::from_path(&paths[0]).expect("error reading bam file {path}");
        let chroms: Vec<String> = bam
            .header()
            .target_names()
            .iter()
            .map(|n| unsafe { str::from_utf8_unchecked(n) }.to_string())
            .collect();
        let output_prefix: PathBuf =
            (output_prefix.to_string_lossy().to_string() + "total-").into();

        let stats = Stat::from_counts(&total_counts, bin_size as usize, ci);
        files::write_stats(stats, output_prefix.clone());
        files::write_errors(&total_counts, output_prefix, chroms);
    }
    Ok(())
}
