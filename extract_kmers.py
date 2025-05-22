import argparse
import csv
import subprocess
import os


def detect_format(path):
    with open(path, 'r') as fh:
        first = fh.read(1)
    if first == '@':
        return 'fastq'
    elif first == '>':
        return 'fasta'
    else:
        raise ValueError('Could not determine file format')


def run_seqkit_sliding(input_file, k):
    """Yield lines of 'name\tsequence' from seqkit sliding."""
    cmd_sliding = ['seqkit', 'sliding', '--step', '1', '--window', str(k), input_file]
    p1 = subprocess.Popen(cmd_sliding, stdout=subprocess.PIPE)
    cmd_fx2tab = ['seqkit', 'fx2tab', '-n', '-s']
    p2 = subprocess.Popen(cmd_fx2tab, stdin=p1.stdout, stdout=subprocess.PIPE, text=True)
    p1.stdout.close()
    for line in p2.stdout:
        yield line.rstrip('\n')
    p2.stdout.close()
    p1.wait()
    p2.wait()


def count_kmers(lines, fmt):
    counts = {}
    read_sets = {} if fmt == 'fastq' else None
    for line in lines:
        name, seq = line.split('\t')
        counts[seq] = counts.get(seq, 0) + 1
        if fmt == 'fastq':
            read_name = name.rsplit(':', 1)[0]
            if seq not in read_sets:
                read_sets[seq] = set()
            read_sets[seq].add(read_name)
    read_counts = {seq: len(rset) for seq, rset in read_sets.items()} if fmt == 'fastq' else None
    return counts, read_counts


def write_csv(counts, read_counts, output_file):
    fieldnames = ['kmer', 'count'] + (['reads'] if read_counts is not None else [])
    with open(output_file, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(fieldnames)
        for kmer, cnt in counts.items():
            row = [kmer, cnt]
            if read_counts is not None:
                row.append(read_counts.get(kmer, 0))
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(description='Extract unique k-mers and their counts using seqkit.')
    parser.add_argument('input', help='Input FASTA/FASTQ file')
    parser.add_argument('-k', '--kmer', type=int, default=30, help='k-mer size (default: 30)')
    parser.add_argument('-o', '--output', default='kmers.csv', help='Output CSV file')
    args = parser.parse_args()

    fmt = detect_format(args.input)
    lines = run_seqkit_sliding(args.input, args.kmer)
    counts, read_counts = count_kmers(lines, fmt)
    write_csv(counts, read_counts, args.output)


if __name__ == '__main__':
    main()
