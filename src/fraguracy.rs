use ndarray::prelude::{Array, ArrayBase, ArrayView5};
use ndarray::{Array5, ArrayViewMut5};
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

fn b() {
    let mut m: Array5<u64> = Array::zeros((2, 50, 5, 5, 6));

    start(m.view_mut());
    start(m.view_mut());
}

fn start(mut counts: ArrayViewMut5<u64>) {
    //let mut counts = Array::zeros((6, 2, 50, 5, 5));
    counts[[0, 0, 20, 2, 2]] += 1;
}

struct ReadInfo {
    read_pos: u8,
    base_q: u8,
    map_q: u8,
    context: u8,
    read: u8, // read 1 or read 2.
    count: u64,
}
