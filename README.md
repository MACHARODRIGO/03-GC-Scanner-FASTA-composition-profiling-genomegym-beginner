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
```

## 📂 Reference Sequences
Example mRNA sequences (downloadable from NCBI) that you can use to test the pipeline:

| Name   | mRNA Accession   | Organism              | Biological function                      | RefSeq / Gene (DNA) |
|:-------|:-----------------|:----------------------|:-----------------------------------------|:--------------------------------------------------------------|
| HBB    | NM_000518.5      | *Homo sapiens*        | Oxygen transport in blood                | [FASTA](https://www.ncbi.nlm.nih.gov/nuccore/NM_000518.5?report=fasta) |
| CYCS   | NM_018947.6      | *Homo sapiens*        | Electron transport in mitochondria       | [FASTA](https://www.ncbi.nlm.nih.gov/nuccore/NM_018947.6?report=fasta) |
| H3F3A  | NM_003548.2      | *Homo sapiens*        | DNA packaging (chromatin, histone H3.3)  | [FASTA](https://www.ncbi.nlm.nih.gov/nuccore/NM_003548.2?report=fasta) |
| INS    | NM_000207.3      | *Homo sapiens*        | Glucose regulation hormone               | [FASTA](https://www.ncbi.nlm.nih.gov/nuccore/NM_000207.3?report=fasta) |
| CAB1   | NM_102733.3      | *Arabidopsis thaliana*| Photosynthesis, pigment-binding          | [FASTA](https://www.ncbi.nlm.nih.gov/nuccore/NM_102733.3?report=fasta) |
