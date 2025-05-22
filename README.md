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


## Locating k-mers in reads

The `locate_kmers.py` script searches a FASTA/FASTQ file for k-mers listed in a CSV file. The file must contain a `kmer` column and may optionally include `count` and `reads` columns produced by `extract_kmers.py`.
When these optional columns are present, only k-mers occurring at least twice and found in at least two reads are processed. Otherwise, all listed k-mers are searched.
For each qualifying k-mer a TSV file containing the locations is written to the output directory (default `kmer_locations`).

```
python locate_kmers.py kmers.csv reads.fastq -o locations
```

Each resulting TSV is named `<kmer>.tsv` and contains the raw output from `seqkit locate` for that k-mer.
