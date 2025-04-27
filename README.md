# kmerlib

**kmerlib** is a lightweight Python toolkit for building and querying *k*-mer indices on large DNA data sets (FASTA genomes or FASTQ reads).

```
📂 kmerlib/
 ├── kmerlib.py       • streaming pass‑1 counter
 ├── filter.py        • keep‑set bitmask (≥X copies) builder
 ├── collect.py       • triple collector  (key, read_id, offset)
 ├── index.py         • chunked binary index writer  (.index.bin)
 ├── query.py         • memory‑mapped API  (K30Index)
📂 tests/
 ├── test_filter.py
 ├── test_index.py
 ├── test_kmer.py
 ├── test_query.py           • pytest checks some functionality (more to come)...
```

---

## Key features

| stage                    | purpose                                                                              | key CLI / API                                                |
| ------------------------ | ------------------------------------------------------------------------------------ | ------------------------------------------------------------ |
| **Pass‑1 counter**       | streams FASTA/FASTQ, rolls canonical 2‑bit hashes, dumps per‑bucket temporary counts | `build_pass1(in_path, is_reads, k, bucket_bits, out_prefix)` |
| **Keep‑set filter**      | selects *k*-mers that appear ≥ *min\_total* times; saves a `.keep.u64` bitmask       | `build_keep_set(prefix, min_total=2)`                        |
| **Occurrence collector** | writes `.occ.tmp` triples `(key, read_id, offset)` plus `.seqmap.bin`                | `collect_occurrences(prefix)`                                |
| **Chunked indexer**      | merges triples into a memory‑mappable `.index.bin` (header + tables)                 | `build_index(prefix, chunk_mb=16)`                           |
| **Query API**            | constant‑time look‑ups: `present`, `count`, `positions`, `count_many`                | `K30Index(map_path)`                                         |

All heavy lifting is done with plain NumPy and memory‑view slices—no C extensions required.

---

## Installation

```bash
# clone and install in editable mode
$ git clone https://github.com/you/kmerlib.git
$ cd kmerlib && pip install -e .[dev]
```

*Requires Python ≥3.8 and NumPy.  The tests depend on **`pytest`**.*

---

## Quick‑start

### 1  Count & filter 30‑mers in a genome

```python
from pathlib import Path
from kmerlib import kmerlib as kl

prefix = Path("lambda_k30")  # output prefix
kl.build_pass1("lambda.fa", is_reads=False, k=30, bucket_bits=14, out_prefix=prefix)
kl.build_keep_set(prefix)                      # default min_total=2
kl.collect_occurrences(prefix)
kl.build_index(prefix)
```

### 2  Query the index

```python
from kmerlib.query import K30Index
idx = K30Index("lambda_k30.index.bin")

print(idx.count("ACGT" * 7))          # how many times does this 30‑mer occur?
print(idx.positions("ACGT" * 7)[:5]) # first five (contig, offset) pairs
```

### 3  Run the test suite

```bash
$ pytest -q kmerlib/tests
```

The tests shipped with *kmerlib* exercise every public stage:

- pass‑1 counting on a small FASTA and FASTQ toy data set,
- keep‑set filtering and occurrence collection,
- index construction and memory‑mapped queries (counts & positions),
- round‑trip integrity: querying a k‑mer found in the input always returns the original coordinates.

These unit tests make no external network calls and finish in <10 seconds on a laptop.  Feel free to add your own cases under `tests/`.

---

## File formats

| suffix          | description                                                            |
| --------------- | ---------------------------------------------------------------------- |
| `.bucket_*.u64` | raw pass‑1 64‑bit counters (one file per bucket)                       |
| `.keep.u64`     | bitmask of surviving keys (1 bit / key)                                |
| `.occ.tmp`      | *uint64* triples `(key, read_id, offset)`                              |
| `.seqmap.bin`   | mapping of `read_id → (file_offset, length)`                           |
| `.index.bin`    | final chunked index: header • key array • offset table • id+pos arrays |

All binaries are little‑endian and **memory‑friendly**: they can be `mmap`‑opened and sliced without loading the whole file.

---

## Performance tips

- Tune `bucket_bits` so that each pass‑1 bucket fits into CPU cache.
- Use `chunk_mb` to limit RAM during the merge stage.
- `count_many()` accepts a Python list/NumPy array of hashes for vectorised look‑ups (>10× faster than looping in Python).

---

## Citing

If you use *kmerlib* in a publication, please cite our spacing‑fingerprint preprint (URL forthcoming).

