# <img src="./img/fastqc-rs-ferris.svg" width=100em alt="fastqc-rs logo" /> peptide finder

![Rust](https://github.com/fxwiegand/fastqc-rs/workflows/Rust/badge.svg)

A fast quality control tool and peptide finder for FASTQ files written in rust inspired by [fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). Results are written to `summary.html` as a self containing html report with visualizations for all statistics.

Available statistics are:
- Read length
- Sequence quality score
- Sequence quality per base
- Sequence content per base
- GC content

## Installation

There are multiple ways to install peptide finder:

#### Source

Download the source code and within the root directory of source run

    cargo install

## Usage

```
pf -q path/to/my_sequence.fastq 
```

Arguments: 

| Parameter                 | Default       | Description   |	
| :------------------------ |:-------------:| :-------------|
| -q --fastq 	       |	-           |The path to the FASTQ file to use
