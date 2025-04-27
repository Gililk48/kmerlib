# tests/test_query.py – Milestone 4 query API validation
import numpy as np, subprocess, urllib.request, pytest, random
from pathlib import Path

from kmerlib import build_pass1
from kmerlib.filter import build_keep_set, collect_occurrences
from kmerlib.index import build_index
from kmerlib.query import K30Index

URL = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/005/845/"
    "GCA_000005845.2_ASM584v2/GCA_000005845.2_ASM584v2_genomic.fna.gz"
)

@pytest.fixture(scope="session")
def idx(tmp_path_factory):
    d = tmp_path_factory.mktemp("m4")
    fasta = d / "lambda.fa"
    if not fasta.exists():
        gz = d / "lambda.fa.gz"
        urllib.request.urlretrieve(URL, gz)
        subprocess.run(["gunzip", "-c", gz], stdout=open(fasta, "wb"), check=True)
    p = d / "lm"
    build_pass1(fasta, is_reads=False, k=30, bucket_bits=10, out_prefix=p)
    build_keep_set(p)
    collect_occurrences(p)
    build_index(p, chunk_mb=16)
    return K30Index(p)


def test_present_and_count(idx):
    # pick 100 random keys from index, check counts match internal counts
    keys = idx._keys
    sample_idx = np.random.choice(len(keys), 100, replace=False)
    for i in sample_idx:
        k = int(keys[i])
        assert idx.present(k)
        assert idx.count(k) == idx._offs[i+1]-idx._offs[i]


def test_count_many(idx):
    sample = np.random.choice(idx._keys, 200, replace=False).astype("<u8")
    counts1 = [idx.count(int(k)) for k in sample]
    counts2 = idx.count_many(sample)
    assert np.all(counts1 == counts2)


def test_positions_consistency(idx):
    # pick 20 keys and compare iterator with direct slice
    for k in np.random.choice(idx._keys, 20, replace=False):
        i = idx._loc(int(k))
        lo,hi = idx._offs[i], idx._offs[i+1]
        ids = idx._ids[lo:hi]; pos = idx._pos[lo:hi]
        it = list(idx.positions(int(k)))
        assert len(it) == hi-lo
        for (rid,p),(rid2,p2) in zip(it, zip(ids,pos)):
            assert rid==rid2 and p==p2
