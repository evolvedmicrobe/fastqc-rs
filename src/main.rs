#[macro_use]
extern crate lazy_static;
mod kmer_processor;
mod process;
use clap::{App, Arg};
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    let matches = App::new("peptide-finder")
        .version("0.1")
        .author("Emily K. Delaney <emilykdelaney@gmail.com> and Felix W. <fxwiegand@wgdnet.de>")
        .about("A quality control tool and combinatorial peptide examination tool for FASTQ files")
        .arg(
            Arg::new("fastq")
                .about("The input FASTQ file to use.")
                .takes_value(true)
                .required(true)
                .short('q')
                .long("fastq"),
        )
        .get_matches();
    crate::process::process(matches.value_of("fastq").unwrap())
}
