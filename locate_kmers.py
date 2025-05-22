import argparse
import csv
import os
import subprocess
import shutil
from tqdm import tqdm
import time


def filter_kmers(kmer_file):
    """Yield kmers to locate.

    If ``count`` and ``reads`` columns are present they are used to filter for
    kmers occurring at least twice in at least two reads. Otherwise all kmers
    in the file are yielded.
    """
    with open(kmer_file, newline='') as fh:
        reader = csv.DictReader(fh)

        # Determine if count/read based filtering can be applied
        use_filters = reader.fieldnames and {
            'count', 'reads'
        }.issubset(reader.fieldnames)

        # Count total lines for progress bar
        fh.seek(0)
        next(fh)  # Skip header
        total_lines = sum(1 for _ in fh)
        fh.seek(0)
        next(fh)  # Skip header again

        for row in tqdm(reader, total=total_lines, desc="Reading k-mers"):
            kmer = row.get('kmer')
            if not kmer:
                continue

            if not use_filters:
                yield kmer
                continue

            try:
                count = int(row.get('count', 0))
                reads = int(row.get('reads', 0))
            except (ValueError, TypeError):
                # skip malformed rows when filtering
                continue

            if count >= 2 and reads >= 2:
                yield kmer


def run_seqkit_locate(kmer, reads_file):
    """Run seqkit locate for a single k-mer."""
    cmd = ['seqkit', 'locate', '-p', kmer, reads_file]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip())
    return proc.stdout


def save_locations(kmer, locations, out_dir):
    """Save k-mer locations to a file."""
    path = os.path.join(out_dir, f"{kmer}.tsv")
    with open(path, 'w') as fh:
        fh.write(locations)


def main():
    parser = argparse.ArgumentParser(description='Locate kmers in reads using seqkit.')
    parser.add_argument(
        'kmers',
        help='CSV file containing a "kmer" column and optional "count" and "reads" columns'
    )
    parser.add_argument('reads', help='Input FASTA/FASTQ file containing reads')
    parser.add_argument('-o', '--outdir', default='kmer_locations', help='Output directory')
    args = parser.parse_args()

    if not shutil.which('seqkit'):
        raise RuntimeError('seqkit not found in PATH')

    os.makedirs(args.outdir, exist_ok=True)

    # Get list of kmers to process
    print("Reading k-mers from file...")
    kmers = list(filter_kmers(args.kmers))
    total_kmers = len(kmers)
    print(f"Found {total_kmers} k-mers to process")

    # Process each k-mer with progress bar
    start_time = time.time()
    for kmer in tqdm(kmers, desc="Locating k-mers"):
        try:
            locs = run_seqkit_locate(kmer, args.reads)
            if locs.strip():  # Only save if we found matches
                save_locations(kmer, locs, args.outdir)
        except Exception as e:
            print(f"\nError processing k-mer {kmer}: {str(e)}")
            continue

    end_time = time.time()
    print(f"\nProcessing completed in {end_time - start_time:.2f} seconds")
    print(f"Average time per k-mer: {(end_time - start_time) / total_kmers:.2f} seconds")


if __name__ == '__main__':
    main()
