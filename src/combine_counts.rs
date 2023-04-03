use crate::fraguracy;
use std::io;
use std::io::BufRead;
use std::path::PathBuf;
use std::string::String;

#[derive(Hash, Debug, PartialOrd, PartialEq, Ord, Eq, Clone)]
pub(crate) struct Count {
    read12: u8,
    orientation: u8,
    read_pos: u8,
    bq_bin: u8,
    context: [char; 2],
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
        self.errors += o.errors;
        self.total += o.total;
    }
}

impl Count {
    fn from_line(s: &str) -> Count {
        let mut sp = s.trim().split('\t');
        eprintln!("splitting: {}", s);
        Count {
            read12: sp.next().expect("not enough columns in line {s}")[1..]
                .parse::<u8>()
                .expect("error parsing int")
                - 1,
            orientation: if sp
                .next()
                .expect("not enough columns in line {s}")
                .chars()
                .nth(1)
                == Some('f')
            {
                0
            } else {
                1
            },
            bq_bin: fraguracy::REVERSE_Q_LOOKUP[sp.next().expect("not enough columns in line {s}")],
            read_pos: sp
                .next()
                .expect("not enough columns in line {s}")
                .parse::<u8>()
                .expect("error parsing int"),
            context: {
                let mut ctx = sp
                    .next()
                    .expect("error getting ctx file from line {s}")
                    .chars();
                [
                    ctx.next().expect("expecting two characters for context"),
                    ctx.next().expect("expecting two characters for context"),
                ]
            },
            total: sp
                .next()
                .expect("not enough columns in line {s}")
                .parse::<u32>()
                .expect("error parsing int"),
            errors: sp
                .next()
                .expect("not enough columns in line {s}")
                .trim()
                .parse::<u32>()
                .expect("error parsing int {sp:?}"),
        }
    }
}

pub(crate) fn combine_counts_main(
    counts_files: Vec<PathBuf>,
    output_path: String,
) -> io::Result<()> {
    let mut counts: std::collections::HashSet<Count> = std::collections::HashSet::new();
    for count_file in counts_files.iter() {
        // open each file and read each line.
        let file = std::fs::File::open(count_file)?;
        let reader = std::io::BufReader::new(file);
        for line in reader.lines() {
            let line = line?;
            let mut c = Count::from_line(&line);
            let entry = counts.get(&c);
            if let Some(entry) = entry {
                c.total += entry.total;
                c.errors += entry.errors;
            }
            counts.insert(c);
        }
    }

    let mut counts: Vec<Count> = counts.into_iter().collect();
    counts.sort();
    for c in counts.iter() {
        dbg!(c);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_line() {
        let line = "r1\tf\t05-19\t0\tAC\t61502\t609";

        let c = Count::from_line(line);
        assert_eq!(c.read12, 0);
        assert_eq!(c.orientation, 1);
        assert_eq!(c.bq_bin, 1);
        assert_eq!(c.read_pos, 0);
        assert_eq!(c.context, ['A', 'C']);
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
            total: 32,
            errors: 1,
        };
        let mut b = a.clone();
        b.errors = 3;
        a += b;

        assert_eq!(a.errors, 4);
    }
}
