use crate::fraguracy::{InnerCounts, Stat};
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

pub(crate) fn write_errors(counts: &InnerCounts, output_prefix: PathBuf, chroms: Vec<String>) {
    let mut errfh = std::fs::File::create(
        output_prefix
            .to_str()
            .expect("error getting output prefix")
            .to_owned()
            + "errors.bed",
    )
    .expect("error opening file!");
    writeln!(errfh, "#chrom\tstart\tend\tbq_bin\tcount").expect("error writing to file");

    for pos in counts.error_positions.keys().sorted() {
        //for (pos, cnt) in (&counts.error_positions).iter() {
        let cnt = counts.error_positions[pos];
        let chrom = &chroms[pos.tid as usize];
        let position = pos.pos;
        let end = position + 1;
        let bqs = crate::fraguracy::Q_LOOKUP[pos.bq_bin as usize];
        writeln!(errfh, "{chrom}\t{position}\t{end}\t{bqs}\t{cnt}")
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
