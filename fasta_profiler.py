import os
import re
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm


# Paths
data_dir = "data"
results_dir = "results"
os.makedirs(results_dir, exist_ok=True)

multifasta = os.path.join(results_dir, "combined_sequences.fasta")
stats_csv = os.path.join(results_dir, "sequence_stats.csv")
meta_csv = os.path.join(results_dir, "metadata.csv")

def combine_fastas(input_folder, output_file):
    """Combine individual FASTA files into a single multifasta."""
    records = []
    for fname in os.listdir(input_folder):
        if fname.endswith(".fasta"):
            records += list(SeqIO.parse(os.path.join(input_folder, fname), "fasta"))
    SeqIO.write(records, output_file, "fasta")

def parse_header_fields(record):
    """
    From a NCBI-like FASTA description, extract:
      - Name: content inside parentheses, e.g., (CAB1)
      - MoleculeType: trailing token after last comma, e.g., 'mRNA'
      - Clean Description: without leading ID, parentheses, or trailing molecule type
    """
    desc = record.description or ""
    # 1) Remove leading ID if present at the beginning
    if desc.startswith(record.id):
        desc = desc[len(record.id):].strip()

    # 2) Extract Name (parentheses)
    name_match = re.search(r"\(([^)]+)\)", desc)
    name = name_match.group(1).strip() if name_match else ""

    # 3) Extract MoleculeType: trailing ", TOKEN" (letters, numbers, dashes, underscores, slashes)
    mol_match = re.search(r",\s*([A-Za-z0-9_\-\/]+)\s*$", desc)
    molecule = mol_match.group(1).strip() if mol_match else ""

    # 4) Build clean description: remove parentheses part and the trailing molecule part
    clean = re.sub(r"\s*\([^)]+\)", "", desc)            # drop "(XYZ)"
    clean = re.sub(r",\s*[A-Za-z0-9_\-\/]+\s*$", "", clean)  # drop ", TOKEN" at end
    clean = clean.strip()

    return name, molecule, clean

def analyze_sequences(fasta_file):
    """Compute composition and GC% for each sequence and return (stats_df, meta_df)."""
    stats_rows, meta_rows = [], []

    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        length = len(seq)
        a, t, c, g = seq.count("A"), seq.count("T"), seq.count("C"), seq.count("G")
        gc_pct = round(100 * (g + c) / length, 2) if length > 0 else 0

        name, molecule, clean_desc = parse_header_fields(record)

        # Stats row (incluye MoleculeType para facilitar análisis)
        stats_rows.append({
            "ID": record.id,
            "Name": name,
            "Length_bp": length,
            "Count_A": a,
            "Count_T": t,
            "Count_C": c,
            "Count_G": g,
            "GC_percent": gc_pct
        })

        # Metadata row (diccionario de datos visible/legible)
        meta_rows.append({
            "ID": record.id,
            "Name": name,
            "MoleculeType": molecule,
            "Description": clean_desc
        })

    return pd.DataFrame(stats_rows), pd.DataFrame(meta_rows)

def plot_summary(stats_df, results_dir):
    """Create a vertical combined figure with sequence length and GC% per gene, with consistent colors and legend."""
    # Paleta de colores (tab10 tiene hasta 10 distintos)
    colors = cm.tab10.colors
    gene_colors = {gene: colors[i % len(colors)] for i, gene in enumerate(stats_df["Name"])}

    fig, axes = plt.subplots(2, 1, figsize=(8, 8))

    # --- Panel A: Length per gene ---
    bars = axes[0].bar(stats_df["Name"],
                       stats_df["Length_bp"],
                       color=[gene_colors[name] for name in stats_df["Name"]],
                       edgecolor="black")
    axes[0].set_xlabel("Gene")
    axes[0].set_ylabel("Length (bp)")
    axes[0].set_title("Sequence length per gene")
    axes[0].tick_params(axis="x", rotation=45)

    for bar in bars:
        height = bar.get_height()
        axes[0].text(bar.get_x() + bar.get_width()/2, height + 20,
                     f"{height}", ha="center", va="bottom", fontsize=8)

    # --- Panel B: GC% per gene ---
    bars2 = axes[1].bar(stats_df["Name"],
                        stats_df["GC_percent"],
                        color=[gene_colors[name] for name in stats_df["Name"]],
                        edgecolor="black")
    axes[1].set_xlabel("Gene")
    axes[1].set_ylabel("GC%")
    axes[1].set_title("GC content per gene")
    axes[1].tick_params(axis="x", rotation=45)

    for bar in bars2:
        height = bar.get_height()
        axes[1].text(bar.get_x() + bar.get_width()/2, height + 1,
                     f"{height:.2f}", ha="center", va="bottom", fontsize=8)


    outfile = os.path.join(results_dir, "summary_plots.png")
    plt.savefig(outfile, dpi=300)
    plt.close()
    print(f"Combined plot saved to {outfile}")

def main():
    # 1) Combinar multifasta
    combine_fastas(data_dir, multifasta)

    # 2) Analizar secuencias
    stats_df, meta_df = analyze_sequences(multifasta)

    # 3) Guardar outputs
    stats_df.to_csv(stats_csv, index=False)
    meta_df.to_csv(meta_csv, index=False)

    print(f"Saved stats to {stats_csv}")
    print(f"Saved metadata to {meta_csv}")

    # 4) Generar gráficos combinados
    plot_summary(stats_df, results_dir)

    # ✅ Preview en consola
    print("\n--- Stats Preview ---")
    print(stats_df.head())
    print("\n--- Metadata Preview ---")
    print(meta_df.head())

if __name__ == "__main__":
    main() # This ensures the script runs when executed directly
# If imported, the main() function won't run automatically.

