use crate::fraguracy;
use std::io;
use std::io::{BufRead, Write};
use std::path::PathBuf;
use std::string::String;

#[derive(Hash, Debug, PartialOrd, PartialEq, Ord, Eq, Clone)]
pub(crate) struct Count {
    read12: u8,
    orientation: u8,
    read_pos: u32,
    bq_bin: u8,
    context: [char; 2],
    homopolymer_dist: i8,
    total: u32,
    errors: u32,
}

impl std::ops::AddAssign<Count> for Count {
    fn add_assign(&mut self, o: Count) {
        assert!(self.read12 == o.read12);
        assert!(self.orientation == o.orientation);
        assert!(self.read_pos == o.read_pos);
        assert!(self.bq_bin == o.bq_bin);
        assert!(self.context == o.context);
        assert!(self.homopolymer_dist == o.homopolymer_dist);
        self.errors += o.errors;
        self.total += o.total;
    }
}

impl Count {
    fn from_line(s: &str, file_name: &str) -> Count {
        let mut sp = s.trim().split('\t');
        Count {
            read12: sp.next().unwrap_or_else(|| panic!("not enough columns in line: {s} from file: {file_name}"))[1..]
                .parse::<u8>()
                .map(|val| val - 1)
                .unwrap_or_else(|e| panic!("error parsing read12 from line: {s} in file: {file_name}, error: {e}")),
            orientation: match sp.next() {
                Some("f") => 0,
                Some("r") => 1,
                _ => panic!("error parsing orientation, expected f or r, in line: {s} from file: {file_name}"),
            },
            bq_bin: fraguracy::REVERSE_Q_LOOKUP[sp.next().unwrap_or_else(|| panic!("not enough columns for bq_bin in line: {s} from file: {file_name}"))],
            read_pos: sp
                .next()
                .unwrap_or_else(|| panic!("not enough columns for read_pos in line: {s} from file: {file_name}"))
                .parse::<u32>()
                .unwrap_or_else(|e| panic!("error parsing read_pos from line: {s} in file: {file_name}, error: {e}")),
            context: {
                let ctx_str = sp
                    .next()
                    .unwrap_or_else(|| panic!("error getting context string from line: {s} from file: {file_name}"));
                let mut ctx_chars = ctx_str.chars();
                [
                    ctx_chars.next().unwrap_or_else(|| panic!("expecting two characters for context, got: {ctx_str} in line: {s} from file: {file_name}")),
                    ctx_chars.next().unwrap_or_else(|| panic!("expecting two characters for context, got: {ctx_str} in line: {s} from file: {file_name}")),
                ]
            },
            homopolymer_dist: sp
                .next()
                .unwrap_or_else(|| panic!("not enough columns for homopolymer_dist in line: {s} from file: {file_name}"))
                .trim()
                .parse::<i8>()
                .unwrap_or_else(|e| panic!("error parsing homopolymer_dist from line: {s} in file: {file_name}, error: {e}")),
            total: sp
                .next()
                .unwrap_or_else(|| panic!("not enough columns for total in line: {s} from file: {file_name}"))
                .parse::<u32>()
                .unwrap_or_else(|e| panic!("error parsing total from line: {s} in file: {file_name}, error: {e}")),
            errors: sp
                .next()
                .unwrap_or_else(|| panic!("not enough columns for errors in line: {s} from file: {file_name}"))
                .trim()
                .parse::<u32>()
                .unwrap_or_else(|e| panic!("error parsing errors from line: {s} in file: {file_name}, error: {e}")),
        }
    }
}

pub(crate) fn combine_counts_main(
    counts_files: Vec<PathBuf>,
    output_path: String,
) -> io::Result<()> {
    let mut counts: std::collections::HashSet<Count> = std::collections::HashSet::new();
    let mut header: String = String::new();
    for count_file in counts_files.iter() {
        // open each file and read each line.
        let file = std::fs::File::open(count_file)?;
        let reader = std::io::BufReader::new(file);
        for (i, line) in reader.lines().enumerate() {
            let line = line?;
            if i == 0 {
                assert!(
                    line.starts_with("read12"),
                    "expecting header line from counts file"
                );
                assert!(
                    line.contains("hp_dist"),
                    "expecting hp_dist in header please run with newer version of fraguracy"
                );
                // take the first 8 columns as the header
                header = line.split('\t').take(8).collect::<Vec<_>>().join("\t");
                continue;
            }
            let mut c = Count::from_line(&line, count_file.to_str().unwrap());
            let entry = counts.take(&c);
            if let Some(entry) = entry {
                c.total += entry.total;
                c.errors += entry.errors;
            }
            counts.insert(c);
        }
    }

    let mut out = std::fs::File::create(output_path)?;
    writeln!(out, "{}", header)?;

    let mut counts: Vec<Count> = counts.into_iter().collect();
    counts.sort();
    for c in counts.iter() {
        writeln!(
            out,
            "r{}\t{}\t{}\t{}\t{}{}\t{}\t{}\t{}",
            c.read12 + 1,
            ['f', 'r'][c.orientation as usize],
            fraguracy::Q_LOOKUP[c.bq_bin as usize],
            c.read_pos,
            c.context[0],
            c.context[1],
            c.homopolymer_dist,
            c.total,
            c.errors
        )?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_line() {
        let line = "r1	f	05-19	0	AC	-1	61502	609";

        let c = Count::from_line(line, "test_file.txt");
        assert_eq!(c.read12, 0);
        assert_eq!(c.orientation, 0);
        assert_eq!(c.bq_bin, 1);
        assert_eq!(c.read_pos, 0);
        assert_eq!(c.context, ['A', 'C']);
        assert_eq!(c.homopolymer_dist, -1);
        assert_eq!(c.total, 61502);
        assert_eq!(c.errors, 609);
    }

    #[test]
    fn test_add_count() {
        let mut a = Count {
            read12: 0,
            orientation: 1,
            bq_bin: 2,
            read_pos: 3,
            context: ['A', 'T'],
            homopolymer_dist: -1,
            total: 32,
            errors: 1,
        };
        let mut b = a.clone();
        b.errors = 3;
        a += b;
        assert_eq!(a.homopolymer_dist, -1);

        assert_eq!(a.errors, 4);
    }
}
