use crate::fraguracy::{InnerCounts, Stat, CONTEXT_TO_CONTEXT2};
use std::string::String;

use itertools::Itertools;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use std::io::Write;

pub(crate) fn write_stats(stats: Vec<Stat>, output_prefix: PathBuf) {
    let header = Stat::header();

    let mut fh = std::fs::File::create(
        output_prefix
            .to_str()
            .expect("error getting output prefix")
            .to_owned()
            + "counts.txt",
    )
    .expect("error opening file!");

    writeln!(fh, "{header}").expect("error writing to file");
    stats
        .iter()
        .for_each(|s| writeln!(fh, "{s}").expect("error writing to file"));
}

pub(crate) fn format_context_counts(counts: [u32; 6]) -> (u32, String) {
    let mut total: u32 = 0;
    let contexts: String = counts
        .iter()
        .enumerate()
        .filter(|(_, &count)| count > 0)
        .map(|(idx, &count)| {
            let context = CONTEXT_TO_CONTEXT2[idx];
            let a = context[0];
            let b = context[1];
            total += count;
            format!("{a}{b}:{count}")
        })
        .collect::<Vec<_>>()
        .join(",");

    (total, contexts)
}

pub(crate) fn write_errors(counts: &InnerCounts, output_prefix: PathBuf, chroms: Vec<String>) {
    let mut errfh = std::fs::File::create(
        output_prefix
            .to_str()
            .expect("error getting output prefix")
            .to_owned()
            + "errors.bed",
    )
    .expect("error opening file!");
    writeln!(errfh, "#chrom\tstart\tend\tbq_bin\tcount\tcontexts").expect("error writing to file");

    for pos in counts.error_positions.keys().sorted() {
        //for (pos, cnt) in (&counts.error_positions).iter() {
        let cnt = counts.error_positions[pos];
        let (total, contexts) = format_context_counts(cnt);
        let chrom = &chroms[pos.tid as usize];
        let position = pos.pos;
        let end = position + 1;
        let bqs = crate::fraguracy::Q_LOOKUP[pos.bq_bin as usize];
        writeln!(
            errfh,
            "{chrom}\t{position}\t{end}\t{bqs}\t{total}\t{contexts}"
        )
        .expect("error writing to error file");
    }
    write_indel_errors(counts, output_prefix, chroms);
}
fn write_indel_errors(counts: &InnerCounts, output_prefix: PathBuf, chroms: Vec<String>) {
    let mut errfh = std::fs::File::create(
        output_prefix
            .to_str()
            .expect("error getting output prefix")
            .to_owned()
            + "indel-errors.bed",
    )
    .expect("error opening indel file!");
    writeln!(errfh, "#chrom\tstart\tend\tcount").expect("error writing to file");
    for pos in counts.indel_error_positions.keys().sorted() {
        let cnt = counts.indel_error_positions[pos];
        let chrom = &chroms[pos.tid as usize];
        let position = pos.pos;
        let end = position + 1;
        writeln!(errfh, "{chrom}\t{position}\t{end}\t{cnt}")
            .expect("error writing to indel-error file");
    }
}

use flate2::read::GzDecoder;
use std::fs::File;

/// Open a file path that may be gzipped.
pub(crate) fn open_file(path: Option<PathBuf>) -> Option<Box<dyn BufRead>> {
    let file = File::open(path.unwrap());
    if file.is_err() {
        eprintln!("error opening file: {}", file.unwrap_err());
        return None;
    }
    let file = file.unwrap();
    let mut buf_file = BufReader::new(file);

    let b = buf_file.fill_buf().expect("error reading from file");
    let gzipped = &b[0..2] == b"\x1f\x8b";

    let reader: Box<dyn BufRead> = if gzipped {
        Box::new(BufReader::new(GzDecoder::new(buf_file)))
    } else {
        Box::new(buf_file)
    };
    Some(reader)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_context_counts() {
        // Test case 1: All counts are non-zero
        let counts1 = [1, 2, 3, 4, 5, 6];
        let (total1, contexts1) = format_context_counts(counts1);
        assert_eq!(total1, 21);
        assert_eq!(contexts1, "AC:1,AG:2,AT:3,CA:4,CG:5,CT:6");

        // Test case 2: Some counts are zero
        let counts2 = [0, 2, 0, 4, 0, 6];
        let (total2, contexts2) = format_context_counts(counts2);
        assert_eq!(total2, 12);
        assert_eq!(contexts2, "AG:2,CA:4,CT:6");

        // Test case 3: All counts are zero
        let counts3 = [0, 0, 0, 0, 0, 0];
        let (total3, contexts3) = format_context_counts(counts3);
        assert_eq!(total3, 0);
        assert_eq!(contexts3, "");

        // Test case 4: Only one non-zero count
        let counts4 = [0, 0, 0, 0, 5, 0];
        let (total4, contexts4) = format_context_counts(counts4);
        assert_eq!(total4, 5);
        assert_eq!(contexts4, "CG:5");
    }
}
