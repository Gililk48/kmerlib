# kmerlib

**kmerlib** is a lightweight Python toolkit for building and querying *k*-mer indices on large DNA data sets (FASTA genomes or FASTQ reads).

Overview
index.py is a Python script designed to process FASTQ (reads) or FASTA (genome) files by extracting k-mers, building an index of their positions, and saving occurrence statistics to a CSV file. The script filters k-mers that appear at least twice (for FASTA) or at least twice in at least two reads (for FASTQ), builds an index mapping k-mers to their positions, and saves the index for later querying. This script is useful for analyzing k-mer distributions and locations in sequencing data, such as in microbiome or genomic studies.
Key Features

Extracts k-mers of a specified length (default: 30) from FASTQ or FASTA files.
Filters k-mers based on occurrence criteria.
Builds an index mapping k-mers to their positions (read/sequence ID and offset).
Saves the index and associated IDs to disk for later querying.
Saves k-mer occurrence statistics to a CSV file.

Requirements
To run index.py, you need the following Python packages installed:

numpy
tqdm

You can install these dependencies using pip:
pip install numpy tqdm

Additionally, ensure you have Python 3.6 or later installed.
Usage
Command-Line Arguments
index.py accepts the following command-line arguments:

input_file (required): Path to the input FASTQ or FASTA file.
--type (required): Type of input file. Choices are:
reads: For FASTQ files (e.g., sequencing reads).
genome: For FASTA files (e.g., genome sequences).


Example usage:

python3 index.py <input_file> --type <reads|genome>

Examples
Example 1: Process a FASTQ File (Reads)
To process a FASTQ file containing sequencing reads:
python3 index.py /path/to/reads.fastq --type reads

This will:

Filter k-mers that appear at least twice in at least two reads.
Build an index mapping k-mers to their positions in the reads.
Save the filtered k-mers to filtered_kmers.txt.
Save the index to index_fastq.pkl and read IDs to ids_fastq.pkl.
Save k-mer stats to filtered_kmer_stats.csv.

Example 2: Process a FASTA File (Genome)
To process a FASTA file containing a genome sequence:
python3 index.py /path/to/genome.fasta --type genome

This will:

Filter k-mers that appear at least twice in the concatenated sequence.
Build an index mapping k-mers to their positions in the sequences.
Save the filtered k-mers to filtered_kmers.txt.
Save the index to index_fasta.pkl and sequence IDs to ids_fasta.pkl.
Save k-mer stats to filtered_kmer_stats.csv.

Output Files
The script generates the following output files:

filtered_kmers.txt: A text file listing all filtered k-mers (one per line).
index_{fastq|fasta}.pkl: A pickled file containing the index (a dictionary mapping k-mers to arrays of positions).
ids_{fastq|fasta}.pkl: A pickled file containing the list of read or sequence IDs corresponding to the indices in the index.
filtered_kmer_stats.csv: A CSV file with k-mer occurrence statistics, with columns:
kmer: The k-mer sequence.
count: Total number of occurrences of the k-mer.
read_id: ID of the read or sequence where the k-mer occurs.
offset: Position (base pair offset) of the k-mer in the read/sequence.



Querying the Saved Index
The saved index (index_{fastq|fasta}.pkl) and IDs (ids_{fastq|fasta}.pkl) can be queried using the query_kmer_index.py script to find the abundance (number of occurrences) and locations of a specific k-mer.
Example: Querying a K-mer
First, run index.py to generate the index:
python3 index.py /path/to/reads.fastq --type reads

Then, use query_kmer_index.py to query a k-mer:
python3 query_kmer_index.py ACGTACGTACGTACGTACGTACGTACGTACGT --type reads

This will output the k-mer’s abundance and locations (read IDs and offsets) based on the saved index.
Notes

K-mer Length: The default k-mer length is 30. To change this, you’d need to modify the k variable in the script.
File Format: The index and IDs are saved using pickle, which is specific to Python. Ensure compatibility with the same Python version when loading the files.
Performance: For large FASTQ or FASTA files, the script may take significant time and memory. The tqdm progress bar provides feedback on the processing status.

License
This script is provided as-is for academic and research purposes. There is no warranty or guarantee of functionality. Use at your own risk.

