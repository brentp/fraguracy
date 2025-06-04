use crate::fraguracy;
use core::cmp::Reverse;
use itertools::Itertools;
use rust_htslib::bgzf;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, Write};
use std::ops::Add;
use std::path::PathBuf;
use std::string::String;

#[derive(Eq, Debug, Default, Clone)]
struct Interval {
    tid: i32,
    chrom: String,
    start: u32,
    end: u32,
    group: u8,
    length: i32,
    count: [u32; 7],
    file_i: u32,
}
struct IntervalHeap {
    // min heap
    h: BinaryHeap<Reverse<Interval>>,
    files: Vec<Box<dyn BufRead>>,
    chom_to_tid: HashMap<String, i32>,
    is_indel: Vec<bool>,
}

fn read_fai(path: PathBuf) -> HashMap<String, i32> {
    let f = File::open(&path);
    let mut h = HashMap::new();
    if let Ok(fu) = f {
        for l in BufReader::new(fu).lines() {
            let l = l.expect("error parsing faidx");
            let chrom = l
                .split('\t')
                .next()
                .expect("expected at least one value per line in faidx");
            if chrom.starts_with('>') {
                log::warn!(
                    "expecting fai, NOT fasta for argument found chrom of {}",
                    chrom
                );
            }
            h.insert(String::from(chrom), h.len() as i32);
        }
    } else {
        panic!("couldn't open file: {:?}", path.to_string_lossy());
    }
    h
}

impl Add<&Interval> for &Interval {
    type Output = Interval;
    fn add(self, other: &Interval) -> Self::Output {
        assert_eq!(self.chrom, other.chrom);
        assert_eq!(self.start, other.start);
        assert_eq!(self.end, other.end);
        assert_eq!(self.group, other.group);
        assert_eq!(self.length, other.length, "indel lengths must be equal");
        let counts = self
            .count
            .iter()
            .zip(other.count.iter())
            .map(|(a, b)| a + b);
        // convert counts to [u32, 7]
        let counts: [u32; 7] = counts
            .collect::<Vec<_>>()
            .try_into()
            .expect("error converting counts");

        Interval {
            chrom: self.chrom.clone(),
            count: counts,
            ..*self
        }
    }
}

impl IntervalHeap {
    fn all_indels(&self) -> bool {
        self.is_indel.iter().all(|&x| x)
    }
    fn new(paths: Vec<PathBuf>, fai_path: PathBuf) -> IntervalHeap {
        let fhs: Vec<Box<dyn BufRead>> = paths
            .iter()
            .map(|p| crate::files::open_file(Some(p.clone())).expect("error opening file"))
            .collect();

        let mut ih = IntervalHeap {
            h: BinaryHeap::new(),
            files: fhs,
            chom_to_tid: read_fai(fai_path),
            is_indel: paths
                .iter()
                .map(|p| p.to_string_lossy().ends_with("indel-errors.bed.gz"))
                .collect(),
        };

        assert!(
            ih.is_indel.iter().all(|&x| x) || ih.is_indel.iter().all(|&x| !x),
            "all files must be either indel error files or base error files, not a mix"
        );

        ih.files
            .iter_mut()
            .enumerate()
            .for_each(|(file_i, fh)| loop {
                // loop to skip '#' comment lines
                let is_indel = ih.is_indel[file_i];
                let mut buf = String::new();
                let line = fh.read_line(&mut buf);
                if line.is_ok() && !buf.starts_with('#') {
                    let r = parse_bed_line(&buf, file_i as u32, &(ih.chom_to_tid), is_indel);
                    if r.is_err() && buf == "" {
                        break;
                    } else {
                        ih.h.push(Reverse(r.unwrap_or_else(|_| {
                            if buf != "" {}
                            panic!("Error parsing first line from file: '{buf}'")
                        })));
                        break;
                    }
                }
            });
        ih
    }
}
fn parse_bed_line(
    line: &str,
    file_i: u32,
    chrom_to_tid: &HashMap<String, i32>,
    is_indel: bool,
) -> Result<Interval, Box<dyn Error>> {
    let toks: Vec<&str> = line.trim().split('\t').collect();
    // can be 6 if combine-errors was already run once.
    let mut iv = if !is_indel {
        let mut iv = Interval {
            tid: 0,
            chrom: String::from(toks[0]),
            start: str::parse::<u32>(toks[1])?,
            end: str::parse::<u32>(toks[2])?,
            group: (*fraguracy::REVERSE_Q_LOOKUP
                .get(toks[3].trim())
                .unwrap_or_else(|| panic!("unknown bq bin: {}", toks[3]))),
            length: 0,
            count: [0; 7],

            file_i,
        };
        // toks[4] is the total count, which we don't need. because we can sum the count from toks[5]

        // parse the counts which appear as, e.g.,
        // AC:1,AG:2,AT:3,CG:4,CT:5,GT:6,NN:0
        // and increment the appropriate index using CONTEXT_LOOKUP from fraguracy.rs
        for s in toks[5].split(',') {
            let (context, count) = s.split(':').collect_tuple().unwrap();
            if context.len() != 2 {
                return Err(
                    format!("expecting two characters for context, found {}", context).into(),
                );
            }
            let mut context = context.chars();
            let a = context.next().unwrap();
            let b = context.next().unwrap();
            let idx = fraguracy::CONTEXT_LOOKUP[&(a as u8, b as u8)];
            iv.count[idx] += count.parse::<u32>().unwrap();
        }
        iv
    } else {
        // indel errors
        let mut iv = Interval {
            tid: 0,
            chrom: String::from(toks[0]),
            start: str::parse::<u32>(toks[1])?,
            end: str::parse::<u32>(toks[2])?,
            group: *fraguracy::REVERSE_Q_LOOKUP
                .get(toks[5].trim())
                .unwrap_or_else(|| panic!("unknown bq bin: {}", toks[5])),
            count: [0; 7],
            length: str::parse::<i32>(toks[4])?,
            file_i,
        };
        // store the count in the first position for indels.
        iv.count[0] = str::parse::<u32>(toks[3])?;
        iv
    };
    iv.tid = *chrom_to_tid
        .get(&iv.chrom)
        .unwrap_or_else(|| panic!("chromosome '{}' not found in fai file", iv.chrom));
    Ok(iv)
}

impl Iterator for IntervalHeap {
    type Item = Interval;

    /// pop an item out and then read in another interval from that file-handle
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(pop_iv) = self.h.pop() {
            let pop_iv = pop_iv.0;
            let file_i = pop_iv.file_i;
            let is_indel = self.is_indel[file_i as usize];
            let fh = &mut self.files[file_i as usize];
            let mut buf = String::new();
            let line_len = &fh.read_line(&mut buf);
            if line_len.is_ok() && *(line_len).as_ref().unwrap() > 0 {
                let r = parse_bed_line(&buf, file_i, &self.chom_to_tid, is_indel);
                if let Ok(iv) = r {
                    self.h.push(Reverse(iv));
                } else {
                    panic!("{:?} line_len: {:?}", r.err().unwrap(), line_len);
                }
            }
            Some(pop_iv)
        } else {
            None
        }
    }
}

impl PartialEq for Interval {
    fn eq(&self, b: &Interval) -> bool {
        self.chrom == b.chrom
            && self.start == b.start
            && self.end == b.end
            && self.group == b.group
            && self.length == b.length
    }
}

impl PartialOrd for Interval {
    #[allow(clippy::non_canonical_partial_ord_impl)]
    fn partial_cmp(&self, b: &Interval) -> Option<Ordering> {
        if self.tid != b.tid {
            return if self.tid < b.tid {
                Some(Ordering::Less)
            } else {
                Some(Ordering::Greater)
            };
        }
        if self.start != b.start {
            return if self.start < b.start {
                Some(Ordering::Less)
            } else {
                Some(Ordering::Greater)
            };
        }

        if self.end != b.end {
            return if self.end < b.end {
                Some(Ordering::Less)
            } else {
                Some(Ordering::Greater)
            };
        }

        if self.length != b.length {
            return if self.length < b.length {
                Some(Ordering::Less)
            } else {
                Some(Ordering::Greater)
            };
        }

        Some(self.group.cmp(&b.group))
    }
}

impl Ord for Interval {
    fn cmp(&self, b: &Interval) -> std::cmp::Ordering {
        self.partial_cmp(b).expect("cmp: not expecting None")
    }
}

pub(crate) fn combine_errors_main(
    paths: Vec<PathBuf>,
    fai_path: PathBuf,
    output_path: String,
) -> io::Result<()> {
    let ih = IntervalHeap::new(paths, fai_path);
    if ih.all_indels() {
        log::info!("all indels");
    }

    // Append .gz if not already present
    let mut output_path = if !output_path.ends_with(".gz") {
        output_path + ".gz"
    } else {
        output_path
    };
    let all_indels = ih.all_indels();

    if all_indels && !output_path.ends_with("indel-errors.bed.gz") {
        log::warn!("all indels, but output path does not end with 'indel-errors.bed.gz'. renaming");
        output_path = output_path.replace(".bed.gz", "indel-errors.bed.gz");
    }

    let mut writer =
        bgzf::Writer::from_path(&output_path).expect("error creating bgzip output file");

    if all_indels {
        writer.write_all(b"#chrom\tstart\tend\t\tcount\tlength\tbq_bin\tn_samples\n")?;
    } else {
        writer.write_all(b"#chrom\tstart\tend\tbq_bin\tcount\tcontexts\tn_samples\n")?;
    }

    for (_, ivs) in &ih
        .into_iter()
        .group_by(|iv| (iv.tid, iv.start, iv.end, iv.group, iv.length))
    {
        let ivs: Vec<Interval> = ivs.into_iter().collect();
        let n = ivs
            .iter()
            .filter(|iv| iv.count.iter().any(|&c| c > 0))
            .count();
        let iv0 = ivs[0].clone();
        let iv = ivs.iter().skip(1).fold(iv0, |acc, iv| &acc + iv);

        let line = if all_indels {
            format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                iv.chrom,
                iv.start,
                iv.end,
                iv.count[0],
                iv.length,
                fraguracy::Q_LOOKUP[iv.group as usize],
                n
            )
        } else {
            let (total_count, context_str) = crate::files::format_context_counts(iv.count);
            format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                iv.chrom,
                iv.start,
                iv.end,
                if iv.group == u8::MAX {
                    "NA"
                } else {
                    fraguracy::Q_LOOKUP[iv.group as usize]
                },
                total_count,
                context_str,
                n
            )
        };
        writer.write_all(line.as_bytes())?;
    }
    log::info!("wrote {}", output_path);
    writer.flush()?;
    Ok(())
}
