use crate::fraguracy;
use core::cmp::Reverse;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::string::String;

#[derive(Eq, Debug)]
struct Interval {
    // NOTE: BSP. might need tid from fai to do sorting if data is too sparse.
    tid: i32,
    chrom: String,
    start: u32,
    end: u32,
    group: u8,
    file_i: u32,
}
struct IntervalHeap {
    // min heap
    h: BinaryHeap<Reverse<Interval>>,
    files: Vec<Box<dyn BufRead>>,
    chom_to_tid: HashMap<String, i32>,
}

fn read_fai(path: PathBuf) -> HashMap<String, i32> {
    let f = File::open(path);
    let mut h = HashMap::new();
    if let Ok(fu) = f {
        for l in BufReader::new(fu).lines() {
            let l = l.expect("error parsing faidx");
            let chrom = l
                .split('\t')
                .into_iter()
                .next()
                .expect("expected at least one value per line in faidx");
            h.insert(String::from(chrom), h.len() as i32);
        }
    }
    h
}

impl IntervalHeap {
    fn new(paths: Vec<PathBuf>, fai_path: PathBuf) -> IntervalHeap {
        let fhs: Vec<Box<dyn BufRead>> = paths
            .iter()
            .map(|p| crate::files::open_file(Some(p.clone())).expect("error opening file"))
            .collect();

        let mut ih = IntervalHeap {
            h: BinaryHeap::new(),
            files: fhs,
            chom_to_tid: read_fai(fai_path),
        };

        ih.files.iter_mut().enumerate().for_each(|(file_i, fh)| {
            let mut buf = String::new();
            let line = fh.read_line(&mut buf);
            if line.is_ok() {
                let r = parse_bed_line(&buf, file_i as u32, &(ih.chom_to_tid));
                ih.h.push(Reverse(r.expect("Error reading first interval from file")))
            }
        });

        ih
    }
}
fn parse_bed_line(
    line: &String,
    file_i: u32,
    chrom_to_tid: &HashMap<String, i32>,
) -> Result<Interval, Box<dyn Error>> {
    let toks: Vec<&str> = line.split('\t').collect();
    let mut iv = Interval {
        tid: 0,
        chrom: String::from(toks[0]),
        start: str::parse::<u32>(toks[1])?,
        end: str::parse::<u32>(toks[2])?,
        group: (*fraguracy::REVERSE_Q_LOOKUP
            .get(toks[3].trim())
            .expect(&format!("unknown bq bin: {}", toks[3]))),
        file_i: file_i,
    };
    iv.tid = chrom_to_tid[&iv.chrom];
    Ok(iv)
}

//fn parse_bed_line(line: &String) -> Result<Interval, std::io::Error> {

impl Iterator for IntervalHeap {
    type Item = Interval;

    /// pop an item out and then read in another interval from that file-handle
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(pop_iv) = self.h.pop() {
            let pop_iv = pop_iv.0;
            let file_i = pop_iv.file_i;
            let fh = &mut self.files[file_i as usize];
            let mut buf = String::new();
            let line = &fh.read_line(&mut buf);
            if line.is_ok() && *(line).as_ref().unwrap() > 0 {
                let r = parse_bed_line(&buf, file_i, &self.chom_to_tid);
                if let Ok(iv) = r {
                    self.h.push(Reverse(iv));
                } else {
                    panic!("{:?}", r.err().unwrap());
                }
            }
            return Some(pop_iv);
        } else {
            return None;
        }
    }
}

impl PartialEq for Interval {
    fn eq(&self, b: &Interval) -> bool {
        return self.chrom == b.chrom
            && self.start == b.start
            && self.end == b.end
            && self.group == b.group;
    }
}

impl PartialOrd for Interval {
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

        Some(self.group.cmp(&b.group))
    }
}

impl Ord for Interval {
    fn cmp(&self, b: &Interval) -> std::cmp::Ordering {
        return self.partial_cmp(b).expect("cmp: not expecting None");
    }
}

pub(crate) fn combine_errors_main(
    paths: Vec<PathBuf>,
    fai_path: PathBuf,
    output_path: String,
) -> io::Result<()> {
    let ih = IntervalHeap::new(paths, fai_path);

    for iv in ih.into_iter() {
        eprintln!("{:?}", iv);
    }
    Ok(())
}
