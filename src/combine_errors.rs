use crate::fraguracy;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::error::Error;
use std::io;
use std::io::BufRead;
use std::path::PathBuf;
use std::string::String;

#[derive(Eq, Ord, Debug)]
struct Interval {
    chrom: String,
    start: u32,
    end: u32,
    group: u8,
    file_i: u32,
}
struct IntervalHeap {
    h: BinaryHeap<Interval>,
    files: Vec<Box<dyn BufRead>>,
}

impl IntervalHeap {
    fn new(paths: Vec<PathBuf>) -> IntervalHeap {
        let fhs: Vec<Box<dyn BufRead>> = paths
            .iter()
            .map(|p| crate::files::open_file(Some(p.clone())).expect("error opening file"))
            .collect();
        let mut ih = IntervalHeap {
            h: BinaryHeap::new(),
            files: fhs,
        };

        ih.files.iter_mut().enumerate().for_each(|(file_i, fh)| {
            let mut buf = String::new();
            let line = &fh.read_line(&mut buf);
            if line.is_ok() {
                let r = parse_bed_line(&buf, file_i as u32);
                ih.h.push(r.expect("Error reading first interval from file"))
            }
        });

        ih
    }
}

//fn parse_bed_line(line: &String) -> Result<Interval, std::io::Error> {
fn parse_bed_line(line: &String, file_i: u32) -> Result<Interval, Box<dyn Error>> {
    let toks: Vec<&str> = line.split('\t').collect();
    let iv = Interval {
        chrom: String::from(toks[0]),
        start: str::parse::<u32>(toks[1])?,
        end: str::parse::<u32>(toks[2])?,
        group: (*fraguracy::REVERSE_Q_LOOKUP
            .get(toks[3])
            .expect(&format!("unknown bq bin: {}", toks[3]))),
        file_i: file_i,
    };
    Ok(iv)
}

impl Iterator for IntervalHeap {
    type Item = Interval;

    /// pop an item out and then read in another interval from that file-handle
    fn next(&mut self) -> Option<Self::Item> {
        let n = self.h.pop();
        if let Some(iv) = self.h.pop() {
            let file_i = iv.file_i;
            let fh = &mut self.files[file_i as usize];
            let mut buf = String::new();
            let line = &fh.read_line(&mut buf);
            if line.is_ok() {
                let r = parse_bed_line(&buf, file_i);
                if let Ok(iv) = r {
                    self.h.push(iv);
                } else {
                    panic!("{:?}", r.err().unwrap());
                }
            }
        }
        return n;
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
        if self.chrom != b.chrom {
            return if self.chrom < b.chrom {
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

pub(crate) fn combine_errors_main(paths: Vec<PathBuf>, output_prefix: String) -> io::Result<()> {
    let ih = IntervalHeap::new(paths);

    for iv in ih.into_iter() {
        eprintln!("{:?}", iv);
    }
    Ok(())
}
