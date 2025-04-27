# tests/test_index.py – Milestone 3 index validation (fixed header handling)
import numpy as np, subprocess, urllib.request, pytest, random, struct
from pathlib import Path

from kmerlib import build_pass1
from kmerlib.filter import build_keep_set, collect_occurrences
from kmerlib.index import build_index

URL = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/005/845/"
    "GCA_000005845.2_ASM584v2/GCA_000005845.2_ASM584v2_genomic.fna.gz"
)

HEADER_FMT = "<8sIIQQ"  # magic, ver, reserved, n_keys, n_occ
HEADER_LEN = struct.calcsize(HEADER_FMT)

@pytest.fixture(scope="session")
def prefix(tmp_path_factory):
    d = tmp_path_factory.mktemp("m3")
    fasta = d / "lambda.fa"
    if not fasta.exists():
        gz = d / "lambda.fa.gz"
        urllib.request.urlretrieve(URL, gz)
        subprocess.run(["gunzip", "-c", gz], stdout=open(fasta, "wb"), check=True)
    p = d / "lm"
    build_pass1(fasta, is_reads=False, k=30, bucket_bits=10, out_prefix=p)
    build_keep_set(p)
    collect_occurrences(p)
    build_index(p, chunk_mb=16)  # small chunk for CI speed
    return p


def _read_header(fp):
    hdr = fp.read(HEADER_LEN)
    magic, ver, res, n_keys, n_occ = struct.unpack(HEADER_FMT, hdr)
    assert magic.startswith(b"KIDX") and ver == 1
    return n_keys, n_occ


def test_offsets_match_counts(prefix):
    idx = prefix.with_suffix(".index.bin")
    with open(idx, "rb") as fp:
        n_keys, n_occ = _read_header(fp)
        keys = np.fromfile(fp, "<u8", n_keys)
        offs = np.fromfile(fp, "<u8", n_keys + 1)
    counts = np.diff(offs)

    count_tbl = np.fromfile(prefix.with_suffix(".count.bin"), "<u8").reshape(-1, 3)
    ref = {int(k): int(t) for k, t, _ in count_tbl}
    for k, c in zip(keys, counts):
        assert ref[k] == c


def test_positions_recoverable(prefix):
    idx = prefix.with_suffix(".index.bin")
    with open(idx, "rb") as fp:
        n_keys, n_occ = _read_header(fp)
        keys = np.fromfile(fp, "<u8", n_keys)
        offs = np.fromfile(fp, "<u8", n_keys + 1)
        ids = np.fromfile(fp, "<u4", n_occ)
        pos = np.fromfile(fp, "<u4", n_occ)

    # sample 30 random keys and ensure slices line up
    sample_ids = random.sample(range(int(n_keys)), min(30, int(n_keys)))
    for i in sample_ids:
        lo, hi = offs[i], offs[i + 1]
        assert hi > lo
        assert hi - lo == len(ids[lo:hi]) == len(pos[lo:hi])
