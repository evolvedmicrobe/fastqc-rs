use memchr::memchr;
use memchr::memmem::{self, Finder};
use needletail::Sequence;
use protein_translate::translate;
use rustc_hash::FxHashMap as HashMap;
use std::convert::TryInto;
use std::fs::File;
use std::io::prelude::*;

// Left flanking AVI-Tag Sequence
// Translates to GLNDIFEAQKIEWHEGGS
const AVI_SEQ: &[u8; 54] = b"GGTCTTAATGATATTTTTGAAGCTCAGAAGATTGAATGGCATGAAGGTGGTAGT";
const AVI_SUFFIX_START: usize = 40;
const AVI_SUFFIX_SEARCH_START: usize = AVI_SUFFIX_START - 5;
lazy_static! {
    static ref PREFIX_FINDER: Finder<'static> = memmem::Finder::new(&AVI_SEQ[40..]);
}
//
#[derive(Default)]
pub(crate) struct CombinatorialPeptideSequenceProcessor {
    pub(crate) cnt_total_reads: usize,
    pub(crate) cnt_valid_front: usize,
    pub(crate) cnt_valid_back: usize,
    pub(crate) peptide_cnts: HashMap<[u8; 7], usize>,
    pub(crate) invalid_cnts: HashMap<[u8; 10], usize>,
}

impl CombinatorialPeptideSequenceProcessor {
    pub(crate) fn new() -> Self {
        Default::default()
    }
    pub(crate) fn process<'a>(&mut self, record: &'a (dyn Sequence<'a> + 'a)) {
        self.cnt_total_reads += 1;
        let full_seq = record.sequence();
        let seq = &full_seq[AVI_SUFFIX_SEARCH_START..full_seq.len()];
        if let Some(i) = PREFIX_FINDER.find(seq) {
            self.cnt_valid_front += 1;
            let remainder = i + PREFIX_FINDER.needle().len();
            let combo_seq = &seq[remainder..];
            let peptide = translate(combo_seq).into_bytes();
            if peptide.len() > 10 {
                if peptide[7] == b'G' && peptide[8] == b'G' && peptide[9] == b'S' {
                    self.cnt_valid_back += 1;
                    let combo_seq: [u8; 7] = peptide[..7].try_into().unwrap();
                    *(self.peptide_cnts.entry(combo_seq).or_insert(0)) += 1;
                } else {
                    let invalid_seq: [u8; 10] = peptide[..10].try_into().unwrap();
                    *(self.invalid_cnts.entry(invalid_seq).or_insert(0)) += 1;
                }
            }
        }
    }

    pub(crate) fn output_results(&mut self) {
        let mut seq_and_cnts: Vec<([u8; 7], usize)> = self.peptide_cnts.drain().collect();
        seq_and_cnts.sort_by_key(|a| std::cmp::Reverse(a.1));
        let mut all_seqs = File::create("All_Results.csv").unwrap();
        let mut no_q_seqs = File::create("NoQ.csv").unwrap();
        let mut outfiles: Vec<File> = (0..7)
            .map(|x| File::create(format!("Q_at_{}.csv", x + 1)).unwrap())
            .collect();
        for (seq, cnt) in seq_and_cnts.into_iter() {
            let formatted =
                format!("{},{}\n", std::str::from_utf8(&seq).unwrap(), cnt).into_bytes();
            all_seqs.write_all(&formatted).unwrap();
            match memchr(b'Q', &seq) {
                Some(i) => {
                    outfiles[i].write_all(&formatted).unwrap();
                }
                None => {
                    no_q_seqs.write_all(&formatted).unwrap();
                }
            }
        }

        // Now output the sequence that didn't parse correctly
        let mut invalid_seq_and_cnts: Vec<([u8; 10], usize)> = self.invalid_cnts.drain().collect();
        invalid_seq_and_cnts.sort_by_key(|a| std::cmp::Reverse(a.1));
        let mut invalid_seqs = File::create("InvalidSeqs.csv").unwrap();
        for (seq, cnt) in invalid_seq_and_cnts.into_iter() {
            let formatted =
                format!("{},{}\n", std::str::from_utf8(&seq).unwrap(), cnt).into_bytes();
            invalid_seqs.write_all(&formatted).unwrap();
        }
    }
}
