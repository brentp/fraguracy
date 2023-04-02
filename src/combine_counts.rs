use crate::fraguracy;
use std::io;
use std::path::PathBuf;
use std::string::String;

#[derive(Hash, Debug, PartialOrd, PartialEq, Ord, Eq)]
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

pub(crate) fn combine_counts_main(counts: Vec<PathBuf>, output_path: String) -> io::Result<()> {
    Ok(())
}

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
}
