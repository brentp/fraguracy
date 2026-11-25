mod combine_counts;
mod combine_errors;
mod depth;
mod homopolymer;
mod pileup;

mod files;
mod fraguracy;
mod lua;

//mod plot;
#[macro_use]
extern crate lazy_static;
use clap::{Parser, Subcommand};
use fraguracy::ConfidenceInterval;
use homopolymer::find_homopolymers;
use linear_map::LinearMap;
use regex::Regex;

use rust_lapper::Lapper;

use std::io::BufRead;

use crate::files::Iv;
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

lazy_static! {
    static ref EMPTY_LAPPER: Lapper<u32, u8> = Lapper::new(Vec::new());
}

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
        about = "combine error bed files from extract. file names are used to determine if they are indel errors or base errors."
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
        #[arg(required = true, help = "bam/cram files to analyze")]
        bams: Vec<PathBuf>,
        #[arg(
            short,
            long,
            default_value_t = String::from("fraguracy-"),
            help = "prefix for output files"
        )]
        output_prefix: String,

        #[arg(short = 'C', long, help = "restrict analysis to this chromosome")]
        chromosome: Option<String>,

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
            short = 'l',
            long,
            help = "optional lua expression to filter reads (once per overlapping pair). returns true to skip pair. e.g. 'return read.flags.secondary or read.flags.supplementary'."
        )]
        lua_expression: Option<String>,

        #[arg(
            short = 'L',
            long,
            help = "optional per-base lua expression (evaluated only if -l passes). has access to `read` and `base` (qpos, bq, gpos, distance_from_5prime, distance_from_3prime, base, is_indel, indel_length). returns true to skip base."
        )]
        lua_base_expression: Option<String>,

        #[arg(
            short,
            long,
            default_value_t = 151,
            help = "indicate the maximum read length in the alignment file"
        )]
        max_read_length: u32,
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
            help = "do not calculate denominator. This can shorten runtime but will also skip the homopolymer distance calculation.",
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

fn get_sample_name(hmap: &HashMap<String, Vec<LinearMap<String, String>>>) -> String {
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

/// Compute 0-based exclusive reference end position directly from `raw_cigar()` without allocating.
#[inline(always)]
fn raw_end_pos(r: &bam::Record) -> i64 {
    // Reference-consuming CIGAR ops: M(0), D(2), N(3), =(7), X(8).
    const REF_CONSUMING_MASK: u32 = (1 << 0) | (1 << 2) | (1 << 3) | (1 << 7) | (1 << 8);
    let mut pos = r.pos();
    for &c in r.raw_cigar() {
        if (REF_CONSUMING_MASK & (1u32 << (c & 0xF))) != 0 {
            pos += (c >> 4) as i64;
        }
    }
    pos
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
            chromosome,
            output_prefix,
            regions,
            exclude_regions,
            lua_expression,
            lua_base_expression,
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
            chromosome,
            PathBuf::from(output_prefix),
            regions,
            exclude_regions,
            lua_expression,
            lua_base_expression,
            bin_size as u32,
            max_read_length,
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

fn read_bed(path: Option<PathBuf>) -> Option<HashMap<String, Lapper<u32, u8>>> {
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

    let mut tree: HashMap<String, Lapper<u32, u8>> = HashMap::new();

    for (chrom, ivs) in bed {
        let mut lapper = Lapper::new(ivs);
        lapper.merge_overlaps();
        tree.insert(chrom, lapper);
    }
    Some(tree)
}

fn get_tree<'a>(
    regions: &'a Option<HashMap<String, Lapper<u32, u8>>>,
    chrom: &String,
) -> Option<&'a Lapper<u32, u8>> {
    if let Some(map) = regions {
        Some(map.get(chrom).unwrap_or(&EMPTY_LAPPER))
    } else {
        None
    }
}

#[allow(clippy::too_many_arguments)]
fn process_bam(
    path: PathBuf,
    fasta_path: Option<PathBuf>,
    regions: Option<PathBuf>,
    chromosome: Option<String>,
    exclude_regions: Option<PathBuf>,
    lua_expression: Option<&str>,
    lua_base_expression: Option<&str>,
    bin_size: u32,
    max_read_length: u32,
    min_mapping_quality: u8,
    min_base_qual: u8,
    reference_as_truth: bool,
    output_prefix: PathBuf,
    no_denominator: bool,
    homopolymer_regex: Option<Regex>,
) -> (fraguracy::InnerCounts, Vec<String>, String) {
    let lua_filter = if lua_expression.is_some() || lua_base_expression.is_some() {
        Some(
            lua::LuaFilter::new(lua_expression, lua_base_expression, mlua::Lua::new())
                .expect("error creating lua interpreter"),
        )
    } else {
        None
    };

    let mut bam = IndexedReader::from_path(&path)
        .unwrap_or_else(|_| panic!("error reading bam file {path:?}"));
    bam.set_threads(2).expect("error setting threads");
    let mut map: FxHashMap<Box<[u8]>, bam::Record> = FxHashMap::default();

    let include_regions = read_bed(regions);
    let exclude_regions = read_bed(exclude_regions);

    let mut ibam = IndexedReader::from_path(&path)
        .unwrap_or_else(|_| panic!("bam file {path:?} must be sorted and indexed"));
    ibam.set_threads(2)
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
        bin_size,
        min_base_qual,
        min_mapping_quality,
        lua_filter,
    );
    let hmap = bam::Header::from_template(bam.header()).to_hashmap();
    let sample_name = get_sample_name(&hmap);
    let is_coord_sorted = hmap
        .get("HD")
        .and_then(|hd| hd.first())
        .and_then(|lm| lm.get("SO"))
        .is_some_and(|so| so == "coordinate");
    log::info!("found sample {sample_name}");

    if !no_denominator {
        counts
            .set_depth_writer(
                &(output_prefix.to_string_lossy().to_string()
                    + &sample_name
                    + "-fraguracy-denominator-depth.bed.gz")
                    .to_string(),
            )
            .expect("error setting depth writer. check permissions/existence of output directory.");
    }

    if let Some(chromosome) = chromosome {
        if let Err(e) = bam.fetch(bam::FetchDefinition::String(chromosome.as_bytes())) {
            log::error!("error fetching chromosome {chromosome}: {e}. iterating over all reads.");
        } else {
            log::info!("limiting analysis to chromosome: \"{chromosome}");
        }
    } else if let Err(e) = bam.fetch(bam::FetchDefinition::All) {
        log::error!("error fetching all reads: {e}");
    }

    let mut n_total = 0;
    let mut n_pairs = 0;
    let chroms: Vec<String> = bam
        .header()
        .target_names()
        .iter()
        .map(|n| unsafe { str::from_utf8_unchecked(n) }.to_string())
        .collect();

    let mut chrom_ctx = fraguracy::ChromCtx {
        chrom: &chroms[0],
        chrom_seq: None,
        include_tree: get_tree(&include_regions, &chroms[0]),
        exclude_tree: get_tree(&exclude_regions, &chroms[0]),
        hp_tree: None,
    };

    let mut last_tid: i32 = -1;
    // Reuse a single `bam::Record` across filtered-out reads and second mates; we only allocate
    // a new `bam::Record` when moving the first mate of a pair into `map`.
    let mut rec = bam::Record::new();
    while let Some(res) = bam.read(&mut rec) {
        n_total += 1;
        res.expect("error parsing read");
        if !fraguracy::filter_read(&rec) {
            continue;
        }
        if rec.is_first_in_template() {
            n_pairs += 1;
        }
        if rec.tid() != last_tid {
            if last_tid != -1 {
                log::info!(
                    "processed chromosome: {} pending orphan reads in map: {}",
                    chroms[last_tid as usize],
                    map.len()
                );
                // Flush remaining depth entries under the previous chromosome name before switching.
                counts.handle_depth(&chroms[last_tid as usize], i64::MAX);
                // In coordinate-sorted mode (with `r.tid() == r.mtid()` enforced by `filter_read`),
                // anything still in `map` after leaving a chromosome is a true orphan.
                if is_coord_sorted {
                    map.clear();
                }
            }
            last_tid = rec.tid();
            let cur_chrom = &chroms[last_tid as usize];
            chrom_ctx.chrom = cur_chrom;

            if include_regions.is_some() {
                chrom_ctx.include_tree = get_tree(&include_regions, cur_chrom);
            }
            if exclude_regions.is_some() {
                chrom_ctx.exclude_tree = get_tree(&exclude_regions, cur_chrom);
            }

            if homopolymer_regex.is_some() || reference_as_truth {
                let seq = fasta
                    .as_ref()
                    .unwrap()
                    .fetch_seq(cur_chrom, 0, i64::MAX as usize)
                    .expect("error fetching sequence from fasta.");
                if let Some(ref re) = homopolymer_regex {
                    chrom_ctx.hp_tree = Some(find_homopolymers(&seq, re));
                }
                if reference_as_truth {
                    chrom_ctx.chrom_seq = Some(seq);
                }
            }
        }

        // Look up by `&[u8]` first via `map.remove` so the second mate of a pair never allocates
        // a qname key or a new `Record`. For coordinate-sorted BAMs we also skip storing reads
        // whose mate cannot overlap them (or that fail MAPQ); name-sorted BAMs keep the generic
        // insert-or-remove path.
        let mate_later = is_coord_sorted && rec.mpos() > rec.pos();
        let mate_earlier = is_coord_sorted && rec.mpos() < rec.pos();
        let found = if mate_later {
            None
        } else {
            map.remove(rec.qname())
        };
        let Some(a) = found else {
            let store = !mate_earlier
                && (!mate_later
                    || (rec.mapq() >= min_mapping_quality && raw_end_pos(&rec) > rec.mpos()));
            if store {
                let key: Box<[u8]> = rec.qname().into();
                map.insert(key, std::mem::replace(&mut rec, bam::Record::new()));
            }
            continue;
        };

        if a.mapq() >= min_mapping_quality && rec.mapq() >= min_mapping_quality {
            // Decode `a`'s CIGAR once here so we can skip non-overlapping pairs before Lua
            // and pass it into `counts.increment` without decoding it a second time.
            let a_cigar = a.cigar();
            if a_cigar.end_pos() > rec.pos() && !counts.skip_pair(&a, &rec) {
                counts.increment(&a, &rec, &a_cigar, &chrom_ctx);
            }
        }
    }
    if last_tid != -1 {
        // Flush any remaining depth entries from the final chromosome at EOF.
        counts.handle_depth(&chroms[last_tid as usize], i64::MAX);
    }
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
    chromosome: Option<String>,
    output_prefix: PathBuf,
    regions: Option<PathBuf>,
    exclude_regions: Option<PathBuf>,
    lua_expression: Option<String>,
    lua_base_expression: Option<String>,
    bin_size: u32,
    max_read_length: u32,
    min_mapping_quality: u8,
    ci: ConfidenceInterval,
    reference_as_truth: bool,
    no_denominator: bool,
    homopolymer_regex: String,
) -> std::io::Result<()> {
    //let args: Vec<String> = env::args().collect();
    let min_base_qual = 5u8;

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
    if reference_as_truth && fasta_path.is_none() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "fasta path must be provided when --reference-as-truth (-t) is set",
        ));
    }

    // Validate lua expressions once up front; each worker thread constructs its own Lua VM.
    if lua_expression.is_some() || lua_base_expression.is_some() {
        lua::LuaFilter::new(
            lua_expression.as_deref(),
            lua_base_expression.as_deref(),
            mlua::Lua::new(),
        )
        .expect("error creating lua interpreter");
    }

    let total_counts = paths
        .par_iter()
        .map(|path| {
            let (c, chroms, sample_name) = process_bam(
                path.clone(),
                fasta_path.clone(),
                regions.clone(),
                chromosome.clone(),
                exclude_regions.clone(),
                lua_expression.as_deref(),
                lua_base_expression.as_deref(),
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    // Assuming Iv is defined in crate::files and is accessible.
    use crate::files::Iv;

    #[test]
    fn test_get_tree_found() {
        // Create a non-empty lapper for "chr1"
        let iv = Iv {
            start: 10,
            stop: 20,
            val: 0,
        };
        let lapper_non_empty = Lapper::new(vec![iv]);
        let mut regions_map: HashMap<String, Lapper<u32, u8>> = HashMap::new();
        regions_map.insert("chr1".to_string(), lapper_non_empty);
        let regions = Some(regions_map);

        let tree = get_tree(&regions, &"chr1".to_string()).unwrap();
        // Query a point that should overlap the interval [10,20]
        let mut iter = tree.find(15, 16);
        assert!(
            iter.next().is_some(),
            "Expected to find an interval for chr1"
        );
    }

    #[test]
    fn test_get_tree_not_found_in_map() {
        // Create a regions map with a lapper only for "chr1"
        let iv = Iv {
            start: 10,
            stop: 20,
            val: 0,
        };
        let lapper_non_empty = Lapper::new(vec![iv]);
        let mut regions_map: HashMap<String, Lapper<u32, u8>> = HashMap::new();
        regions_map.insert("chr1".to_string(), lapper_non_empty);
        let regions = Some(regions_map);

        // Looking up "chr2" should yield the empty lapper.
        let tree = get_tree(&regions, &"chr2".to_string()).unwrap();
        assert!(
            tree.find(0, 100).next().is_none(),
            "Expected empty lapper for non-existent chromosome"
        );
    }

    #[test]
    fn test_get_tree_no_regions() {
        // When regions is None, we expect the function to return None.
        let regions: Option<HashMap<String, Lapper<u32, u8>>> = None;
        let tree = get_tree(&regions, &"any".to_string());
        assert!(tree.is_none(), "Expected None when regions is None");
    }

    #[test]
    fn test_raw_end_pos_matches_cigar() {
        use rust_htslib::bam::record::{Cigar, CigarString};
        let cigar = CigarString(vec![
            Cigar::SoftClip(3),
            Cigar::Match(10),
            Cigar::Ins(2),
            Cigar::Del(4),
            Cigar::Equal(5),
            Cigar::Diff(1),
            Cigar::RefSkip(100),
            Cigar::HardClip(5),
        ]);
        let mut r = bam::Record::new();
        // Query-consuming ops: 3S + 10M + 2I + 5= + 1X = 21 bp
        r.set(b"q1", Some(&cigar), &[b'A'; 21], &[30; 21]);
        r.set_pos(500);
        // Ref-consuming ops: 10M + 4D + 5= + 1X + 100N = 120 bp -> end_pos = 620
        assert_eq!(raw_end_pos(&r), r.cigar().end_pos());
        assert_eq!(raw_end_pos(&r), 620);
    }

    #[test]
    fn test_reference_as_truth_no_denominator_requires_fasta() {
        let err = extract_main(
            vec![],
            None,
            None,
            PathBuf::from("test-"),
            None,
            None,
            None,
            None,
            3,
            151,
            50,
            ConfidenceInterval::AgrestiCoull,
            true,
            true,
            homopolymer::HP_REGEX.to_string(),
        )
        .unwrap_err();
        assert_eq!(err.kind(), std::io::ErrorKind::InvalidInput);
    }
}
