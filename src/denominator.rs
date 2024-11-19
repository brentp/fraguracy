use rust_htslib::bam::{Read, Reader};
use std::io::Result;
use std::path::PathBuf;

pub(crate) fn denominator_main(
    bam_path: PathBuf,
    min_homopolymer_length: u8,
    max_homopolymer_distance: u8,
    fasta: PathBuf,
    min_mapping_quality: u8,
    output_prefix: String,
) -> Result<()> {
    let mut bam = Reader::from_path(bam_path)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
    bam.set_reference(fasta)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
    bam.set_threads(1)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;

    let contigs = bam.header().reference_names();

    Ok(())
}
