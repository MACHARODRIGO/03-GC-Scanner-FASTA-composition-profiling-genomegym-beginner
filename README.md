# 03 – GC Scanner: FASTA Composition Profiling (Beginner)

## 🧬 Goal
Load DNA sequences from a FASTA file, compute nucleotide counts (A/T/C/G) and GC content percentage, then export tidy summary tables to CSV.

---

## 📌 Why it matters
- GC content is a basic descriptor of genomic sequences that varies across organisms and genomic regions.  
- Nucleotide counts provide useful QC metrics and features for downstream machine learning models.

---

## ⚙️ Prerequisites
- Python ≥ 3.9  
- Packages: [Biopython](https://biopython.org/), [pandas](https://pandas.pydata.org/)  

Install dependencies:
```bash
pip install biopython pandas