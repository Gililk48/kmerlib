# tests/test_filter.py – validate Milestone 2 outputs (keep.u64 scheme)
import numpy as np, pytest, random, subprocess, urllib.request
from pathlib import Path

from kmerlib import build_pass1
from kmerlib.filter import build_keep_set, collect_occurrences

URL = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/005/845/"
    "GCA_000005845.2_ASM584v2/GCA_000005845.2_ASM584v2_genomic.fna.gz"
)

@pytest.fixture(scope="session")
def dataset(tmp_path_factory):
    """Run pass-1 and pass-2 on the real lambda-phage genome and return prefix Path."""
    d = tmp_path_factory.mktemp("m2")
    fasta = d / "lambda.fa"
    if not fasta.exists():
        gz = d / "lambda.fa.gz"
        urllib.request.urlretrieve(URL, gz)
        subprocess.run(["gunzip", "-c", gz], stdout=open(fasta, "wb"), check=True)
    prefix = d / "lm"
    build_pass1(fasta, is_reads=False, k=30, bucket_bits=10, out_prefix=prefix)
    build_keep_set(prefix)
    collect_occurrences(prefix)
    return prefix


def test_keep_nonempty_sorted(dataset):
    keep = np.fromfile(dataset.with_suffix(".keep.u64"), dtype="<u8")
    assert len(keep) > 0 and np.all(np.diff(keep) > 0)


def test_occurrence_file_nonempty(dataset):
    occ_path = dataset.with_suffix(".occ.tmp")
    assert occ_path.stat().st_size > 0 and occ_path.stat().st_size % 16 == 0


def test_kept_keys_consistency(dataset):
    # Load keep-set and occurrences
    keep = set(int(x) for x in np.fromfile(dataset.with_suffix(".keep.u64"), dtype="<u8"))
    occ  = np.memmap(dataset.with_suffix(".occ.tmp"), dtype=[("key","<u8"),("id","<u4"),("pos","<u4")])

    # 1. Every key in occurrences must belong to keep-set
    assert set(int(k) for k in np.unique(occ["key"])) <= keep

    # 2. Spot-check 50 kept keys: occurrence count == total in count.bin
    counts = {int(k): int(t) for k,t,_ in np.fromfile(dataset.with_suffix(".count.bin"), dtype="<u8").reshape(-1,3) if k in keep}
    if counts:
        sample = random.sample(list(counts.keys()), min(50, len(counts)))
        for k in sample:
            assert (occ["key"] == k).sum() == counts[k]
