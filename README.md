# ⚛️ GC Scanner: FASTA Composition Profiling

A small, reproducible exercise that:
1) merges multiple FASTA files into a single **multifasta**,
2) computes nucleotide counts and GC content,
3) exports tidy CSV tables, and
4) draws quick QC plots (length & GC%).

---

## 🎯 Learning Objectives  
📥 Load sequence data using Biopython’s `SeqIO`.  
🧮 Count nucleotides (A/T/C/G) and compute GC% content.  
📤 Export results as clean CSV tables.  
📚 Store relevant metadata for future analyses.  
📂 Combine input sequences into a single multifasta file.  
📊 Visualize sequence length distribution to check variability.

## 📌 Why it matters
- GC content is a basic descriptor of genomic sequences that varies across organisms and genomic regions.  
- Nucleotide counts provide useful QC metrics and features for downstream machine learning models.

---

## ⚙️ Prerequisites
- Python ≥ 3.9  
- Packages: [Biopython](https://biopython.org/), [pandas](https://pandas.pydata.org/),   [matplotlib](https://matplotlib.org/stable/)

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


**Note:** These sequences are biologically valid and error-free.  
They vary in length and function — useful for testing your GC scanner and visualizations.

## 🚀 Usage

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/MACHARODRIGO/03-GC-Scanner-FASTA-composition-profiling-genomegym-beginner.git](https://github.com/MACHARODRIGO/03-GC-Scanner-FASTA-composition-profiling-genomegym-beginner.git)
    cd 03-GC-Scanner-FASTA-composition-profiling-genomegym-beginner
    ```
2.  **Download prerequisite data:**
    - Create a directory named `data/`.
    - Download the example FASTA files from the links in the "Reference Sequences" table and save them inside the `data/` folder.
    - **Note:** Ensure the filenames end with `.fasta` (e.g., `HBB.fasta`).
3.  **Run the script:**
    - Use the command-line arguments to specify your input and output directories.
    - The --input flag defaults to data, --output defaults to results, and --pattern defaults to *.fasta.
    - 
    ```bash
     python fasta_profiler.py --input data --output results --pattern "*.fasta"
    ```
5.  **Explore the results:**
    - The script will create a `results/` directory containing:
      - `combined_sequences.fasta`
      - `sequence_stats.csv`
      - `metadata.csv`
      - `summary_plots.png`

## 📊 Example Results

After running the script, the `results/` folder will contain the following files:

### **`sequence_stats.csv` preview**
```csv
ID,Name,Length_bp,Count_A,Count_T,Count_C,Count_G,GC_percent
NM_102733.3,CAB1,1230,316,330,294,290,47.48
NM_018947.6,CYCS,5432,1548,1645,1030,1209,41.22
...
```

### **`metadata.csv` preview**
```csv
ID,Name,MoleculeType,Description
NM_102733.3,CAB1,mRNA,Arabidopsis thaliana chlorophyll A/B binding protein 1
NM_018947.6,CYCS,mRNA,"Homo sapiens cytochrome c, somatic"
...
```

## 💡 Next Steps & Community Challenges

* **Error Handling:** Modify the `analyze_sequences` function to handle sequences with characters other than A, T, C, G (e.g., N for unknown bases). You could print a warning or skip the sequence.
* **Advanced Visualization:** Add a third subplot to `plot_summary` that shows the A/T ratio or a stacked bar chart of all nucleotide counts.
* **Batch Processing:** Adapt the script to handle a large number of FASTA files more efficiently, perhaps by using `os.walk` to traverse a directory tree.
* **Testing:** Write unit tests for the functions, especially `parse_header_fields` and `analyze_sequences`, to ensure they work correctly with different inputs.
