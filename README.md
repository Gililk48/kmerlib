K-mer Indexer (index.py)

Extracts k-mers from FASTQ/FASTA files, builds an index, and saves stats for genomic analysis.
Features

Filters k-mers (length 30) appearing ≥2 times (FASTA) or ≥2 times in ≥2 reads (FASTQ).
Builds and saves an index of k-mer positions.
Saves stats to CSV.

Installation
Requires Python 3.6+ and:
pip install numpy tqdm

Usage
python index.py <input_file> --type <reads|genome>


input_file: Path to FASTQ/FASTA file.
--type: reads (FASTQ) or genome (FASTA).

Output

filtered_kmers.txt: Filtered k-mers.
index_{fastq|fasta}.pkl: Saved index.
ids_{fastq|fasta}.pkl: Saved read/sequence IDs.
filtered_kmer_stats.csv: K-mer stats (kmer, count, read_id, offset).

Querying
Query the index with query_kmer_index.py:
python query_kmer_index.py ACGTACGT... --type reads

License
MIT License. See LICENSE.
