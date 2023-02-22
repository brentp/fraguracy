use rust_htslib::bam::Record;
use std::collections::HashMap;
use std::rc::Rc;

lazy_static! {
    pub(crate) static ref contextLookup: HashMap<(u8, u8), u8> = HashMap::from([
        (('C' as u8, 'A' as u8), 0u8),
        (('G' as u8, 'T' as u8), 0u8),
        (('C' as u8, 'G' as u8), 1u8),
        (('G' as u8, 'C' as u8), 1u8),
        (('C' as u8, 'T' as u8), 2u8),
        (('G' as u8, 'A' as u8), 2u8),
        (('T' as u8, 'A' as u8), 3u8),
        (('A' as u8, 'T' as u8), 3u8),
        (('T' as u8, 'C' as u8), 4u8),
        (('A' as u8, 'G' as u8), 4u8),
        (('T' as u8, 'G' as u8), 5u8),
        (('A' as u8, 'C' as u8), 5u8),
    ]);
}

pub(crate) fn filter_read(r: &Rc<Record>) -> bool {
    r.tid() == r.mtid()
        && r.tid() >= 0
        && !r.is_unmapped()
        && !r.is_mate_unmapped()
        && (r.pos() - r.mpos()).abs() < 1000
        && !r.is_supplementary()
        && !r.is_secondary()
        && !r.is_duplicate()
        && !r.is_quality_check_failed()
}

enum Context {
    CA = 0,
    CG,
    CT,
    TA,
    TC,
    TG,
}
struct ReadInfo {
    readPos: u8,
    baseQ: u8,
    mapQ: u8,
    context: Context,
}

struct Fragment {
    read1: ReadInfo,
    read2: ReadInfo,
}
