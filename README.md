# K‑mer Indexer

**`index.py`** – A lightweight command‑line tool for extracting _k_-mers from FASTQ/FASTA files, building a position index, and exporting simple per‑k‑mer statistics for downstream genomic analysis.

---

## ✨ Features

* Filters 30‑mers that appear **≥ 2** times in a genome (FASTA) or **≥ 2** times in **≥ 2** reads (FASTQ).
* Builds a fast lookup index that maps each retained _k_-mer to `(read|seq_id, offset)` positions.
* Persists the index and associated metadata to disk for rapid re‑use.
* Exports a CSV with summary counts so you can explore the data in Excel/R/Python without re‑processing the raw files.

---

## 🚀 Installation

### Prerequisites

* Python **3.6+**

### Install dependencies

```bash
pip install numpy tqdm
```

---

## 📖 Usage

```bash
python index.py <input_file> --type <reads|genome>
```

| Argument | Description |
|----------|-------------|
| `<input_file>` | Path to a **FASTQ** or **FASTA** file. |
| `--type` | `reads` (for FASTQ) or `genome` (for FASTA). |

### Examples

Index a paired‑end read set in **reads.fastq**:

```bash
python index.py reads.fastq --type reads
```

Index a reference genome in **ref.fa**:

```bash
python index.py ref.fa --type genome
```

---

## 📂 Output

| File | When produced | What it contains |
|------|---------------|------------------|
| `filtered_kmers.txt` | Always | One retained 30‑mer per line. |
| `index_fastq.pkl` / `index_fasta.pkl` | Always | `dict[kmer] -> List[(id, offset)]` pickle. |
| `ids_fastq.pkl` / `ids_fasta.pkl` | Always | List of read / sequence identifiers. |
| `filtered_kmer_stats.csv` | Always | CSV with columns: `kmer,count,id,offset`. |

---

## 🔍 Querying the Index

Retrieve location(s) of a specific 30‑mer after the index has been built:

```bash
python query_kmer_index.py ACGTACGTACGTACGTACGTACGTACGTAC --type reads
```

The script prints matching `(id, offset)` tuples and summary counts to STDOUT.

---

## 🤝 Contributing

Pull requests are welcome! If you find a bug or have a feature request, feel free to open an issue first.

---

## ⚖️ License

This project is licensed under the **MIT License** – see the [LICENSE](LICENSE) file for details.

