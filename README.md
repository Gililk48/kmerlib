# kmerlib

This repository contains a simple script for extracting unique k-mers from
FASTA or FASTQ files using `seqkit`.

## Requirements

- [seqkit](https://bioinf.shenwei.me/seqkit/) must be installed and available
  on your `$PATH`.
- Python 3

## Usage

```
python extract_kmers.py <input.fasta|input.fastq> [-k KMER_SIZE] [-o OUTPUT_CSV]
```

The script will output a CSV file (default `kmers.csv`) containing the k-mer and
its total count. When processing FASTQ files the CSV also includes the number of
reads in which each k-mer appears.

```
# example
python extract_kmers.py reads.fastq -k 30 -o kmers.csv
```

