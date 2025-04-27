# kmerlib/index.py – Milestone 3 (external sort + index writer)
# ---------------------------------------------------------------------------
# Public API
#   build_index(prefix, chunk_mb=1024)
#   chunk_sort(prefix, chunk_bytes)
#   merge_chunks(prefix)
#
# Index.bin layout (little‑endian)
#   magic 8 B  (b"KIDXv1\0\0")
#   ver   4 B  (uint32 = 1)
#   res   4 B  (reserved)
#   nKey  8 B  (#unique k‑mers)
#   nOcc  8 B  (#occurrence rows)
#   key  [nKey] uint64  sorted
#   off  [nKey+1] uint64  offsets into id/pos arrays
#   id   [nOcc] uint32
#   pos  [nOcc] uint32
# ---------------------------------------------------------------------------
from __future__ import annotations

import heapq, json, struct, sys
from pathlib import Path

import numpy as np

REC_DTYPE = np.dtype([("key","<u8"),("id","<u4"),("pos","<u4")])
REC_SIZE  = 16
MAGIC     = b"KIDXv1\0\0"
HEADER_FMT= "<8sIIQQ"  # magic, ver, res, nKey, nOcc

# ───────────────────── phase 1: chunk sort ────────────────────────────────

def chunk_sort(prefix: str | Path, chunk_bytes: int = 1<<30):
    prefix = Path(prefix)
    occ = prefix.with_suffix(".occ.tmp")
    n_total = occ.stat().st_size // REC_SIZE
    chunk_elems = max(1, chunk_bytes // REC_SIZE)
    mmap = np.memmap(occ, dtype=REC_DTYPE)
    n_chunks = 0
    for start in range(0, n_total, chunk_elems):
        end = min(start+chunk_elems, n_total)
        chunk = np.sort(mmap[start:end], order="key")
        out = prefix.with_suffix(f".sorted.{n_chunks:03d}.bin")
        chunk.tofile(out)
        n_chunks += 1
        print(f"[M3] chunk {n_chunks}: {(end-start)*REC_SIZE/1e6:.1f} MB")
    return n_chunks

# ───────────────────── phase 2: k‑way merge ───────────────────────────────

def _stream(path: Path):
    arr = np.memmap(path, dtype=REC_DTYPE)
    return [arr,0]

def merge_chunks(prefix: str | Path):
    prefix = Path(prefix)
    parts = sorted(prefix.parent.glob(prefix.name + ".sorted.*.bin"))
    streams = [_stream(p) for p in parts]
    heap: list[tuple[int,int]] = []
    for i,(a,_) in enumerate(streams):
        if len(a): heap.append((int(a[0]["key"]), i))
    heapq.heapify(heap)

    keys, offs, ids, pos = [], [0], [], []
    cur = None
    while heap:
        k,s = heapq.heappop(heap)
        arr, idx = streams[s]
        rec = arr[idx]
        if k!=cur:
            if cur is not None: offs.append(len(ids))
            keys.append(k); cur=k
        ids.append(int(rec["id"]))
        pos.append(int(rec["pos"]))
        idx+=1; streams[s][1]=idx
        if idx<len(arr): heapq.heappush(heap,(int(arr[idx]["key"]),s))
    offs.append(len(ids))

    idx_path = prefix.with_suffix(".index.bin")
    with open(idx_path,"wb") as fp:
        fp.write(struct.pack(HEADER_FMT, MAGIC,1,0,len(keys),len(ids)))
        np.array(keys,dtype="<u8").tofile(fp)
        np.array(offs,dtype="<u8").tofile(fp)
        np.array(ids ,dtype="<u4").tofile(fp)
        np.array(pos ,dtype="<u4").tofile(fp)
    print(f"[M3] index → {idx_path} ({idx_path.stat().st_size/1e6:.1f} MB)")

    for p in parts: p.unlink()

    meta_path=prefix.with_suffix(".bkidx.json")
    with open(meta_path) as fh: meta=json.load(fh)
    meta.update({"index_file":str(idx_path),"n_keys":len(keys),"n_occ":len(ids)})
    with open(meta_path,"w") as fh: json.dump(meta,fh,indent=2)

# ───────────────────── wrapper & CLI ───────────────────────────────────────

def build_index(prefix: str | Path, chunk_mb:int=1024):
    chunk_sort(prefix, chunk_mb<<20)
    merge_chunks(prefix)

def cli(argv=None):
    import argparse
    ap=argparse.ArgumentParser(description="Milestone 3: build index.bin")
    ap.add_argument("prefix"); ap.add_argument("--chunk-mb",type=int,default=1024)
    args=ap.parse_args(argv)
    build_index(args.prefix,args.chunk_mb)

if __name__=="__main__":
    cli()
