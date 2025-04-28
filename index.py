# index.py: Functions for k-mer filtering, indexing, and stats saving from FASTQ or FASTA files.

from collections import defaultdict, Counter
import numpy as np
from tqdm import tqdm
import time
import csv
import argparse

def count_lines(file_path):
    """
    Count the total number of lines in a file.

    Args:
        file_path (str): Path to the input file.

    Returns:
        int: Total number of lines in the file.
    """
    with open(file_path, 'r') as f:
        return sum(1 for _ in f)

def get_filtered_kmers_fastq(file_path, k=30):
    """
    Identify k-mers from a FASTQ file that appear at least twice in at least two reads.

    Args:
        file_path (str): Path to the FASTQ file.
        k (int): Length of k-mers (default: 30).

    Returns:
        set: Set of k-mers that appear >=2 times in >=2 reads.

    Side Effect:
        Writes filtered k-mers to 'filtered_kmers.txt'.
    """
    start_time = time.time()
    
    # Calculate total reads for progress bar (4 lines per read in FASTQ)
    total_lines = count_lines(file_path)
    total_reads = total_lines // 4

    # Track how many reads each k-mer appears in with >=2 occurrences
    kmer_read_count = defaultdict(int)

    with open(file_path, 'r') as f:
        pbar = tqdm(total=total_reads, desc="FASTQ: filtering k-mers")
        while True:
            try:
                next(f)  # Skip read ID line
                seq = next(f).strip()  # Read sequence
                next(f)  # Skip '+' line
                next(f)  # Skip quality scores

                # Extract k-mers from the sequence
                if len(seq) >= k:
                    kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
                    kmer_freq = Counter(kmers)
                    # Count reads where k-mer appears >=2 times
                    for kmer, freq in kmer_freq.items():
                        if freq >= 2:
                            kmer_read_count[kmer] += 1

                pbar.update(1)
            except StopIteration:
                break
        pbar.close()

    # Filter k-mers that appear in >=2 reads
    filtered_kmers = {kmer for kmer, count in kmer_read_count.items() if count >= 2}
    
    # Save filtered k-mers to file
    with open('filtered_kmers.txt', 'w') as f:
        for kmer in filtered_kmers:
            f.write(f"{kmer}\n")
    
    elapsed_time = time.time() - start_time
    print(f"FASTQ k-mer filtering completed in {elapsed_time:.2f} seconds")
    return filtered_kmers

def get_filtered_kmers_fasta(file_path, k=30):
    """
    Identify k-mers from a FASTA file that appear at least twice in the concatenated sequence.

    Args:
        file_path (str): Path to the FASTA file.
        k (int): Length of k-mers (default: 30).

    Returns:
        set: Set of k-mers that appear >=2 times in the concatenated sequence.

    Side Effect:
        Writes filtered k-mers to 'filtered_kmers.txt'.
    """
    start_time = time.time()
    
    # Read and concatenate all sequences from the FASTA file
    sequences = []
    current_seq = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    sequences.append(''.join(current_seq))
                    current_seq = []
            else:
                current_seq.append(line)
        if current_seq:
            sequences.append(''.join(current_seq))
    
    full_sequence = ''.join(sequences)
    
    # Extract and count k-mers
    if len(full_sequence) >= k:
        kmers = [full_sequence[i:i+k] for i in range(len(full_sequence) - k + 1)]
        kmer_freq = Counter(kmers)
        filtered_kmers = {kmer for kmer, count in kmer_freq.items() if count >= 2}
    else:
        filtered_kmers = set()
    
    # Save filtered k-mers to file
    with open('filtered_kmers.txt', 'w') as f:
        for kmer in filtered_kmers:
            f.write(f"{kmer}\n")
    
    elapsed_time = time.time() - start_time
    print(f"FASTA k-mer filtering completed in {elapsed_time:.2f} seconds")
    return filtered_kmers

def build_index_fastq(file_path, filtered_kmers, k=30):
    """
    Build an index for filtered k-mers from a FASTQ file, mapping k-mers to their positions.

    Args:
        file_path (str): Path to the FASTQ file.
        filtered_kmers (set): Set of k-mers to index.
        k (int): Length of k-mers (default: 30).

    Returns:
        tuple: (index, read_ids)
            - index: Dictionary mapping k-mers to arrays of (read_index, offset).
            - read_ids: List of read IDs corresponding to read indices.
    """
    start_time = time.time()
    
    index = defaultdict(list)
    read_ids = []

    # Calculate total reads for progress bar
    total_lines = count_lines(file_path)
    total_reads = total_lines // 4

    with open(file_path, 'r') as f:
        pbar = tqdm(total=total_reads, desc="FASTQ: building index")
        while True:
            try:
                id_line = next(f).strip()  # Read ID line
                seq = next(f).strip()  # Sequence
                next(f)  # Skip '+' line
                next(f)  # Skip quality scores

                read_id = id_line
                read_index = len(read_ids)
                read_ids.append(read_id)

                # Index k-mers in the sequence
                if len(seq) >= k:
                    for offset in range(len(seq) - k + 1):
                        kmer = seq[offset:offset + k]
                        if kmer in filtered_kmers:
                            index[kmer].append((read_index, offset))

                pbar.update(1)
            except StopIteration:
                break
        pbar.close()

    # Convert lists to numpy arrays for efficient storage
    for kmer in index:
        index[kmer] = np.array(
            index[kmer],
            dtype=[('read_index', 'uint32'), ('offset', 'uint16')]
        )

    elapsed_time = time.time() - start_time
    print(f"FASTQ index building completed in {elapsed_time:.2f} seconds")
    return index, read_ids

def build_index_fasta(file_path, filtered_kmers, k=30):
    """
    Build an index for filtered k-mers from a FASTA file, mapping k-mers to their positions.

    Args:
        file_path (str): Path to the FASTA file.
        filtered_kmers (set): Set of k-mers to index.
        k (int): Length of k-mers (default: 30).

    Returns:
        tuple: (index, seq_ids)
            - index: Dictionary mapping k-mers to arrays of (seq_index, offset).
            - seq_ids: List of sequence IDs corresponding to sequence indices.
    """
    start_time = time.time()
    
    index = defaultdict(list)
    seq_ids = []

    with open(file_path, 'r') as f:
        pbar = tqdm(desc="FASTA: building index", total=count_lines(file_path))
        current_seq = []
        current_id = None
        line_count = 0
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq and current_id:
                    seq = ''.join(current_seq)
                    seq_index = len(seq_ids)
                    seq_ids.append(current_id)
                    if len(seq) >= k:
                        for offset in range(len(seq) - k + 1):
                            kmer = seq[offset:offset + k]
                            if kmer in filtered_kmers:
                                index[kmer].append((seq_index, offset))
                current_id = line
                current_seq = []
            else:
                current_seq.append(line)
            line_count += 1
            pbar.update(1)
        
        if current_seq and current_id:
            seq = ''.join(current_seq)
            seq_index = len(seq_ids)
            seq_ids.append(current_id)
            if len(seq) >= k:
                for offset in range(len(seq) - k + 1):
                    kmer = seq[offset:offset + k]
                    if kmer in filtered_kmers:
                        index[kmer].append((seq_index, offset))
        pbar.close()

    # Convert lists to numpy arrays for efficient storage
    for kmer in index:
        index[kmer] = np.array(
            index[kmer],
            dtype=[('read_index', 'uint32'), ('offset', 'uint16')]
        )

    elapsed_time = time.time() - start_time
    print(f"FASTA index building completed in {elapsed_time:.2f} seconds")
    return index, seq_ids

def query_kmer(index, read_ids, kmer):
    """
    Query the index for a k-mer's occurrences.

    Args:
        index (dict): Index mapping k-mers to arrays of (read_index, offset).
        read_ids (list): List of read or sequence IDs.
        kmer (str): k-mer to query.

    Returns:
        tuple: (count, results)
            - count (int): Number of occurrences of the k-mer.
            - results (list): List of (read_id, offset) tuples for each occurrence.
    """
    if kmer in index:
        occurrences = index[kmer]
        count = len(occurrences)
        results = [(read_ids[occ['read_index']], occ['offset']) for occ in occurrences]
        return count, results
    return 0, []

def save_kmer_stats(index, read_ids, filtered_kmers):
    """
    Query each filtered k-mer and save occurrence stats to CSV.

    Args:
        index (dict): Index mapping k-mers to arrays of (read_index, offset).
        read_ids (list): List of read or sequence IDs.
        filtered_kmers (set): Set of k-mers to query and save.

    Side Effect:
        Writes stats to 'filtered_kmer_stats.csv' with columns: kmer, count, read_id, offset.
    """
    with open('filtered_kmer_stats.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['kmer', 'count', 'read_id', 'offset'])
        
        for kmer in tqdm(filtered_kmers, desc="Querying k-mers and saving to CSV"):
            count, occurrences = query_kmer(index, read_ids, kmer)
            for read_id, offset in occurrences:
                writer.writerow([kmer, count, read_id, offset])

def main():
    """
    Main function to process FASTQ or FASTA files for k-mer indexing.

    Parses command-line arguments, filters k-mers, builds an index, and saves stats.
    """
    parser = argparse.ArgumentParser(description="Process FASTQ (reads) or FASTA (genome) for k-mer indexing.")
    parser.add_argument('input_file', type=str, help="Path to the input FASTQ or FASTA file")
    parser.add_argument('--type', type=str, choices=['reads', 'genome'], required=True,
                        help="Type of input: 'reads' for FASTQ, 'genome' for FASTA")
    args = parser.parse_args()

    file_path = args.input_file
    k = 30

    if args.type == 'reads':
        print("Processing FASTQ file (reads)")
        filtered_kmers = get_filtered_kmers_fastq(file_path, k)
        print(f"Number of filtered k-mers: {len(filtered_kmers)}")
        index, read_ids = build_index_fastq(file_path, filtered_kmers, k)
    else:  # args.type == 'genome'
        print("Processing FASTA file (genome)")
        filtered_kmers = get_filtered_kmers_fasta(file_path, k)
        print(f"Number of filtered k-mers: {len(filtered_kmers)}")
        index, read_ids = build_index_fasta(file_path, filtered_kmers, k)

    save_kmer_stats(index, read_ids, filtered_kmers)

if __name__ == "__main__":
    main()