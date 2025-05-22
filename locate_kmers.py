import argparse
import csv
import os
import subprocess
import shutil


def filter_kmers(kmer_file):
    """Yield kmers that occur at least twice in at least two reads."""
    with open(kmer_file, newline='') as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            kmer = row.get('kmer')
            try:
                count = int(row.get('count', 0))
                reads = int(row.get('reads', 0))
            except ValueError:
                continue
            if count >= 2 and reads >= 2:
                yield kmer


def run_seqkit_locate(kmer, reads_file):
    cmd = ['seqkit', 'locate', '-p', kmer, reads_file]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip())
    return proc.stdout


def save_locations(kmer, locations, out_dir):
    path = os.path.join(out_dir, f"{kmer}.tsv")
    with open(path, 'w') as fh:
        fh.write(locations)


def main():
    parser = argparse.ArgumentParser(description='Locate kmers in reads using seqkit.')
    parser.add_argument('kmers', help='CSV file with kmer, count and reads columns')
    parser.add_argument('reads', help='Input FASTA/FASTQ file containing reads')
    parser.add_argument('-o', '--outdir', default='kmer_locations', help='Output directory')
    args = parser.parse_args()

    if not shutil.which('seqkit'):
        raise RuntimeError('seqkit not found in PATH')

    os.makedirs(args.outdir, exist_ok=True)

    for kmer in filter_kmers(args.kmers):
        locs = run_seqkit_locate(kmer, args.reads)
        save_locations(kmer, locs, args.outdir)


if __name__ == '__main__':
    main()
