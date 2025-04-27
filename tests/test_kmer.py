# tests/test_kmer.py – pytest harness for Pass‑1 counter
# ---------------------------------------------------------------------------
# Usage:  pytest -q tests/test_kmer.py
# ---------------------------------------------------------------------------
# Validates the pass‑1 counter on the *real* lambda‑phage genome. All imports
# now reference the renamed `kmerlib` module.
# ---------------------------------------------------------------------------

import random
import subprocess
import urllib.request
from pathlib import Path

import numpy as np
import pytest

from kmerlib import (
    BUCKET_STRUCT,
    build_pass1,
    parse_fastq,
    parse_fasta,
    rolling_canonical_kmers,
    _DNA2BIT,
    _RC2BIT,
)

# ────────────────────────── data fixture ────────────────────────────────────
LAMBDA_FASTA_URL = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/005/845/"
    "GCA_000005845.2_ASM584v2/GCA_000005845.2_ASM584v2_genomic.fna.gz"
)

@pytest.fixture(scope="session")
def lambda_fasta(tmp_path_factory):
    d = tmp_path_factory.mktemp("lambda")
    gz_path = d / "lambda.fa.gz"
    fasta_path = d / "lambda.fa"
    if not fasta_path.exists():
        urllib.request.urlretrieve(LAMBDA_FASTA_URL, gz_path)
        subprocess.run(["gunzip", "-c", gz_path], check=True, stdout=open(fasta_path, "wb"))
    return fasta_path

# ────────────────────────── build counts fixture ───────────────────────────

@pytest.fixture(scope="session")
def count_file(lambda_fasta, tmp_path_factory):
    out = tmp_path_factory.mktemp("idx") / "lm"
    build_pass1(lambda_fasta, is_reads=False, k=30, bucket_bits=10, out_prefix=out)
    return Path(f"{out}.count.bin"), lambda_fasta

# ────────────────────────── slow Python reference counter ──────────────────

def _count_kmers_python(fasta: Path, k: int = 30):
    from collections import defaultdict
    counts = defaultdict(int)
    mask = (1 << (2 * k)) - 1
    high2 = 2 * (k - 1)
    for _, seq in parse_fasta(fasta):
        fwd = rev = 0; n = 0
        for ch in seq:
            bits = _DNA2BIT[ord(ch)]
            fwd = ((fwd << 2) | bits) & mask
            rev = (rev >> 2) | (_RC2BIT[bits] << high2)
            n += 1
            if n >= k:
                counts[min(fwd, rev)] += 1
    return counts

# ────────────────────────── tests ───────────────────────────────────────────

def test_total_kmers(count_file):
    count_bin, fasta = count_file
    buf = np.fromfile(count_bin, dtype="<u8").reshape(-1, 3)
    total_from_file = int(buf[:, 1].sum())
    seq_len = sum(len(seq) for _, seq in parse_fasta(fasta))
    assert total_from_file == seq_len - 30 + 1


def test_unique_keys(count_file):
    count_bin, _ = count_file
    buf = np.fromfile(count_bin, dtype="<u8").reshape(-1, 3)
    keys = np.sort(buf[:, 0])
    assert np.all(np.diff(keys) != 0)


@pytest.mark.slow
def test_sampled_counts(count_file):
    count_bin, fasta = count_file
    buf = np.fromfile(count_bin, dtype="<u8").reshape(-1, 3)
    sample = random.sample(range(buf.shape[0]), k=min(100, buf.shape[0]))
    py_counts = _count_kmers_python(fasta)
    for idx in sample:
        key, total, _ = buf[idx]
        assert py_counts[key] == total
