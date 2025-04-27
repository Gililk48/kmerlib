# kmerlib.py – Pass‑1 k‑mer counter (Milestone 1)
# ---------------------------------------------------------------------------
# * FASTA / FASTQ streaming parser (plain‑text only)
# * 2‑bit DNA encoder & canonical‑30‑mer rolling hash
# * Radix‑bucketed counter writing `{prefix}.count.bin` + `{prefix}.bkidx.json`
# ---------------------------------------------------------------------------
# This is **library code only**.  Tests live in `tests/test_kmer.py`.
# ---------------------------------------------------------------------------

from __future__ import annotations

import argparse
import json
import os
import struct
import sys
from pathlib import Path
from typing import Iterator, Tuple

from tqdm import tqdm  # progress bars

# ────────────────────────── 2‑bit helpers ──────────────────────────────────
_DNA2BIT = {ord("A"): 0, ord("C"): 1, ord("G"): 2, ord("T"): 3}
_RC2BIT  = {0: 3, 1: 2, 2: 1, 3: 0}
BUCKET_STRUCT = struct.Struct("<QQQ")  # key, total, distinct

__all__ = [
    "BUCKET_STRUCT",
    "build_pass1",
    "parse_fastq",
    "parse_fasta",
    "rolling_canonical_kmers",
]

# ────────────────────────── FASTQ / FASTA parsers ──────────────────────────

def parse_fastq(path: str | Path) -> Iterator[Tuple[str, str]]:
    with open(path, "rt", encoding="ascii", errors="strict") as fh:
        while True:
            hdr = fh.readline()
            if not hdr:
                break
            seq = fh.readline(); fh.readline(); fh.readline()
            yield hdr.rstrip("\n")[1:], seq.rstrip("\n")

def parse_fasta(path: str | Path) -> Iterator[Tuple[str, str]]:
    name, buf = None, []
    with open(path, "rt", encoding="ascii", errors="strict") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:].rstrip("\n"), []
            else:
                buf.append(line.rstrip("\n"))
        if name is not None:
            yield name, "".join(buf)

# ────────────────────────── rolling canonical‑k‑mer generator ──────────────

def rolling_canonical_kmers(seq: str, k: int = 30):
    if len(seq) < k:
        return
    mask  = (1 << (2 * k)) - 1
    high2 = 2 * (k - 1)
    fwd = rev = n = 0
    for ch in seq.encode("ascii"):
        bits = _DNA2BIT.get(ch, 255)
        if bits == 255:                # N or invalid
            n = fwd = rev = 0
            continue
        fwd = ((fwd << 2) | bits) & mask
        rev = (rev >> 2) | (_RC2BIT[bits] << high2)
        n += 1
        if n >= k:
            yield min(fwd, rev)

# ────────────────────────── Pass‑1 builder (lazy‑FD + tqdm) ────────────────

def _estimate_reads(path: Path, is_reads: bool):
    if not is_reads:
        return None
    try:
        # fast line count via os.stat + average line len would be nicer, but fine
        with open(path, "rb") as fh:
            n_lines = sum(1 for _ in fh)
        return n_lines // 4
    except Exception:
        return None


def build_pass1(
    input_path: str | Path,
    is_reads: bool,
    k: int,
    bucket_bits: int,
    out_prefix: str | Path,
    max_entries: int = 5_000_000,
):
    input_path = Path(input_path)
    nbuckets   = 1 << bucket_bits
    tmp_paths  = [Path(f"{out_prefix}.bucket{b:04x}.tmp") for b in range(nbuckets)]
    bucket_maps = [{} for _ in range(nbuckets)]

    # tqdm bar for parsing
    n_est = _estimate_reads(input_path, is_reads)
    pbar_parse = tqdm(total=n_est, unit="read", desc="Pass‑1 parse/hash", dynamic_ncols=True)

    def _append(bid: int, payload: bytes):
        with open(tmp_paths[bid], "ab") as fh:
            fh.write(payload)

    # secondary bar for big flush
    flush_bar = None

    def flush(bid: int):
        nonlocal flush_bar
        mp = bucket_maps[bid]
        if not mp:
            return
        if flush_bar is None:
            flush_bar = tqdm(total=sum(len(m) for m in bucket_maps),
                             desc="Flushing buckets", unit="entry", dynamic_ncols=True)
        for key, (tot, dist) in mp.items():
            _append(bid, BUCKET_STRUCT.pack(key, tot, dist))
        flush_bar.update(len(mp))
        mp.clear()

    seq_iter = parse_fastq if is_reads else parse_fasta

    for ridx, (_, seq) in enumerate(seq_iter(input_path)):
        seen = set()
        for key in rolling_canonical_kmers(seq, k):
            bucket = key >> (60 - bucket_bits)
            m      = bucket_maps[bucket]
            tot, dist = m.get(key, (0, 0))
            tot += 1
            if key not in seen:
                dist += 1; seen.add(key)
            m[key] = (tot, dist)
            if len(m) >= max_entries:
                flush(bucket)
        pbar_parse.update(1)
    pbar_parse.close()

    # final flush
    for b in range(nbuckets):
        flush(b)
    if flush_bar is not None:
        flush_bar.close()

    # concatenate → count.bin (with progress)
    count_path = Path(f"{out_prefix}.count.bin")
    concat_bar = tqdm(tmp_paths, desc="Concat buckets", unit="file", dynamic_ncols=True)
    with open(count_path, "wb") as out:
        for p in concat_bar:
            if p.exists():
                out.write(p.read_bytes()); p.unlink()
    concat_bar.close()

    # metadata
    meta = {
        "k": k,
        "bucket_bits": bucket_bits,
        "nbuckets": nbuckets,
        "input": str(input_path),
        "reads": is_reads,
        "count_file": str(count_path),
    }
    with open(f"{out_prefix}.bkidx.json", "w", encoding="utf-8") as fh:
        json.dump(meta, fh, indent=2)

# ────────────────────────── CLI entrypoint ────────────────────────────────

def _cli(argv=None):
    ap = argparse.ArgumentParser(description="Pass‑1 k‑mer counter")
    ap.add_argument("--input", required=True)
    ap.add_argument("--type", choices=["reads", "genome"], required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--k", type=int, default=30)
    ap.add_argument("--bucket-bits", type=int, default=14)
    args = ap.parse_args(argv)
    build_pass1(args.input, args.type == "reads", args.k, args.bucket_bits, args.out)

if __name__ == "__main__":
    _cli()
