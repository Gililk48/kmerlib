# kmerlib/filter.py – Milestone 2 (circl‑import‑free)
# ---------------------------------------------------------------------------
# Outputs for prefix produced by pass‑1:
#   • prefix.keep.u64   – sorted uint64 array of accepted k‑mer keys
#   • prefix.seqmap.bin – 0‑terminated read/contig names → uint32 IDs
#   • prefix.occ.tmp    – unsorted triples (key:uint64, id:uint32, pos:uint32)
# ---------------------------------------------------------------------------
from __future__ import annotations

import json, struct, sys
from pathlib import Path

import numpy as np

# Import *only* from kmerlib.core via package root to avoid circular import.
from kmerlib import (
    BUCKET_STRUCT,
    parse_fastq,
    parse_fasta,
    rolling_canonical_kmers,
)

OCC_STRUCT = struct.Struct("<QI I")  # 16‑byte occurrence rows

# ────────────────────────── helpers ────────────────────────────────────────

def _meta(prefix: Path):
    with open(prefix.with_suffix(".bkidx.json"), "r") as fh:
        return json.load(fh)


def _iter_counts(path: Path):
    row = BUCKET_STRUCT.size
    with open(path, "rb") as fh:
        while chunk := fh.read(row):
            yield BUCKET_STRUCT.unpack(chunk)  # key, tot, dist

# ────────────────────────── keep‑set builder ───────────────────────────────

def build_keep_set(prefix: str | Path, tot_thr: int = 2, distinct_thr: int | None = None):
    """Create prefix.keep.u64 with accepted k‑mer keys.

    * For read datasets (meta["reads"] == True) use default distinct_thr=2
      → keeps k‑mers seen in at least 2 distinct reads.
    * For genome/contig datasets (meta["reads"] == False) we automatically
      relax to distinct_thr = 1 because there is only one sequence ID.
    """
    prefix = Path(prefix)
    meta   = _meta(prefix)
    if distinct_thr is None:
        distinct_thr = 2 if meta["reads"] else 1

    keys = [k for k,t,d in _iter_counts(Path(meta["count_file"]))
            if t >= tot_thr and d >= distinct_thr]

    keys.sort()
    np.array(keys, dtype="<u8").tofile(prefix.with_suffix(".keep.u64"))
    print(f"[M2] keep‑set → {prefix.with_suffix('.keep.u64')}  ({len(keys)} keys)")

# ────────────────────────── occurrence collection ──────────────────────────

def collect_occurrences(prefix: str | Path):
    prefix = Path(prefix)
    meta   = _meta(prefix)
    k      = meta["k"]
    is_reads = meta["reads"]

    keep_keys = set(int(x) for x in np.fromfile(prefix.with_suffix(".keep.u64"), "<u8"))

    seqmap_f = prefix.with_suffix(".seqmap.bin").open("wb")
    occ_f    = prefix.with_suffix(".occ.tmp").open("wb")

    seq_iter = parse_fastq if is_reads else parse_fasta
    rid = 0
    for idx, (name, seq) in enumerate(seq_iter(meta["input"])):
        seqmap_f.write(name.encode("ascii") + b"\0")
        for pos, key in enumerate(rolling_canonical_kmers(seq, k)):
            if key in keep_keys:
                occ_f.write(OCC_STRUCT.pack(key, rid, pos))
        rid += 1
        if idx and idx % 100_000 == 0:
            print(f"[M2] scanned {idx} records", file=sys.stderr)

    seqmap_f.close(); occ_f.close()
    print(f"[M2] occurrences → {occ_f.name}")
    print(f"[M2] seq‑map     → {seqmap_f.name}")

# ────────────────────────── CLI ────────────────────────────────────────────

def cli(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="Milestone 2: build keep‑set + occurrences")
    ap.add_argument("prefix")
    args = ap.parse_args(argv)
    build_keep_set(args.prefix)
    collect_occurrences(args.prefix)

if __name__ == "__main__":
    cli()
