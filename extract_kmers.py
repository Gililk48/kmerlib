import argparse
import csv
import subprocess
import os
from tqdm import tqdm
from collections import defaultdict
import tempfile


def detect_format(path):
    with open(path, 'r') as fh:
        first = fh.read(1)
    if first == '@':
        return 'fastq'
    elif first == '>':
        return 'fasta'
    else:
        raise ValueError('Could not determine file format')


def process_batch(batch_lines, fmt):
    """Process a batch of k-mer lines."""
    counts = defaultdict(int)
    read_sets = defaultdict(set) if fmt == 'fastq' else None
    
    for line in batch_lines:
        name, seq = line.rstrip('\n').split('\t')
        counts[seq] += 1
        if fmt == 'fastq':
            read_name = name.rsplit(':', 1)[0]
            read_sets[seq].add(read_name)
    
    read_counts = {seq: len(rset) for seq, rset in read_sets.items()} if fmt == 'fastq' else None
    return dict(counts), read_counts


def count_kmers_stream(input_file, k, fmt, batch_size=1000000):
    """Count k-mers as they are being extracted from the input file."""
    print(f"Starting k-mer extraction with k={k} from {input_file}")
    
    # Initialize counters
    total_counts = defaultdict(int)
    total_read_counts = defaultdict(int) if fmt == 'fastq' else None
    
    # Run seqkit sliding
    cmd_sliding = ['seqkit', 'sliding', '--step', '1', '--window', str(k), input_file]
    print(f"Running command: {' '.join(cmd_sliding)}")
    p1 = subprocess.Popen(cmd_sliding, stdout=subprocess.PIPE)
    
    # Run fx2tab to get name and sequence
    cmd_fx2tab = ['seqkit', 'fx2tab', '-n', '-s']
    print(f"Running command: {' '.join(cmd_fx2tab)}")
    p2 = subprocess.Popen(cmd_fx2tab, stdin=p1.stdout, stdout=subprocess.PIPE, text=True)
    p1.stdout.close()
    
    # Get total number of sequences for progress bar
    print("Counting total sequences...")
    cmd_count = ['seqkit', 'seq', '-n', input_file]
    total_lines = int(subprocess.check_output(['wc', '-l'], input=subprocess.check_output(cmd_count)).decode().strip())
    print(f"Found {total_lines} sequences to process")
    
    # Process k-mers in batches
    print("Starting k-mer processing...")
    batch = []
    kmer_count = 0
    batch_count = 0
    
    for line in tqdm(p2.stdout, total=total_lines, desc="Processing k-mers"):
        batch.append(line)
        kmer_count += 1
        
        if len(batch) >= batch_size:
            batch_count += 1
            print(f"\nProcessing batch {batch_count} ({kmer_count} k-mers so far)")
            
            # Process batch
            counts, read_counts = process_batch(batch, fmt)
            
            # Merge results
            for seq, count in counts.items():
                total_counts[seq] += count
            if fmt == 'fastq':
                for seq, count in read_counts.items():
                    total_read_counts[seq] = max(total_read_counts[seq], count)
            
            # Clear batch
            batch = []
            
            # Print memory usage
            print(f"Found {len(total_counts)} unique k-mers so far")
    
    # Process final batch if any
    if batch:
        batch_count += 1
        print(f"\nProcessing final batch {batch_count}")
        counts, read_counts = process_batch(batch, fmt)
        for seq, count in counts.items():
            total_counts[seq] += count
        if fmt == 'fastq':
            for seq, count in read_counts.items():
                total_read_counts[seq] = max(total_read_counts[seq], count)
    
    print("\nFinished processing k-mers")
    p2.stdout.close()
    p1.wait()
    p2.wait()
    
    print(f"Found {len(total_counts)} unique k-mers")
    return dict(total_counts), dict(total_read_counts) if fmt == 'fastq' else None


def write_csv(counts, read_counts, output_file):
    """Write results to CSV."""
    print(f"Writing results to {output_file}...")
    fieldnames = ['kmer', 'count'] + (['reads'] if read_counts is not None else [])
    
    # Write in batches to reduce memory usage
    batch_size = 100000
    items = list(counts.items())
    
    with open(output_file, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(fieldnames)
        
        for i in tqdm(range(0, len(items), batch_size), desc="Writing results"):
            batch = items[i:i + batch_size]
            rows = []
            for kmer, cnt in batch:
                row = [kmer, cnt]
                if read_counts is not None:
                    row.append(read_counts.get(kmer, 0))
                rows.append(row)
            writer.writerows(rows)
    
    print("Finished writing results")


def main():
    parser = argparse.ArgumentParser(description='Extract unique k-mers and their counts using seqkit.')
    parser.add_argument('input', help='Input FASTA/FASTQ file')
    parser.add_argument('-k', '--kmer', type=int, default=30, help='k-mer size (default: 30)')
    parser.add_argument('-o', '--output', default='kmers.csv', help='Output CSV file')
    parser.add_argument('-b', '--batch-size', type=int, default=1000000,
                      help='Number of k-mers to process in each batch (default: 1000000)')
    args = parser.parse_args()

    print(f"Processing file: {args.input}")
    fmt = detect_format(args.input)
    print(f"Detected format: {fmt}")
    
    counts, read_counts = count_kmers_stream(args.input, args.kmer, fmt, args.batch_size)
    write_csv(counts, read_counts, args.output)
    print("Done!")


if __name__ == '__main__':
    main()
