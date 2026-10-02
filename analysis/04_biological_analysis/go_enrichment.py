import pandas as pd
import numpy as np
import gseapy as gp

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt

# 1. LOAD ANNOTATED PYDESEQ2 RESULTS
results = pd.read_csv(
    "results/tables/pydeseq2_annotated_results.csv"
)

print("Genes loaded:")
print(len(results))

# 2. KEEP GENES WITH VALID SYMBOLS
results = results.dropna(
    subset=["gene_symbol"]
).copy()

# 3. DEFINE BACKGROUND GENES
# All annotated genes that were tested by PyDESeq2
background_genes = (
    results["gene_symbol"]
    .dropna()
    .drop_duplicates()
    .tolist()
)

print("\nBackground genes:")
print(len(background_genes))

# 4. DEFINE SIGNIFICANT DEGs
significant = results[
    (results["padj"] < 0.05) &
    (results["log2FoldChange"].abs() >= 1)
].copy()

# 5. SEPARATE UP AND DOWN GENES
up_genes = (
    significant.loc[
        significant["log2FoldChange"] >= 1,
        "gene_symbol"
    ]
    .dropna()
    .drop_duplicates()
    .tolist()
)

down_genes = (
    significant.loc[
        significant["log2FoldChange"] <= -1,
        "gene_symbol"
    ]
    .dropna()
    .drop_duplicates()
    .tolist()
)

print("\nUpregulated genes:")
print(len(up_genes))

print("\nDownregulated genes:")
print(len(down_genes))

# 6. DOWNLOAD GO BIOLOGICAL PROCESS 2026
print("\nDownloading GO Biological Process 2026...")

go_bp = gp.get_library(
    name="GO_Biological_Process_2026",
    organism="Mouse"
)

print("GO terms downloaded:")
print(len(go_bp))

# 7. ENRICHMENT FOR UPREGULATED GENES
print("\nRunning enrichment for upregulated genes...")

up_enr = gp.enrichr(
    gene_list=up_genes,
    gene_sets=go_bp,
    background=background_genes,
    outdir=None,
    cutoff=0.05,
    no_plot=True
)

up_results = up_enr.results.copy()

# 8. ENRICHMENT FOR DOWNREGULATED GENES
print("\nRunning enrichment for downregulated genes...")

down_enr = gp.enrichr(
    gene_list=down_genes,
    gene_sets=go_bp,
    background=background_genes,
    outdir=None,
    cutoff=0.05,
    no_plot=True
)

down_results = down_enr.results.copy()

# 9. SAVE FULL RESULTS
up_results.to_csv(
    "results/tables/go_upregulated.csv",
    index=False
)

down_results.to_csv(
    "results/tables/go_downregulated.csv",
    index=False
)

# 10. FILTER SIGNIFICANT GO TERMS
up_sig = up_results[
    up_results["Adjusted P-value"] < 0.05
].copy()

down_sig = down_results[
    down_results["Adjusted P-value"] < 0.05
].copy()

print("\n================================")
print("GO ENRICHMENT SUMMARY")
print("================================")

print("\nSignificant UP GO terms:")
print(len(up_sig))

print("\nSignificant DOWN GO terms:")
print(len(down_sig))

# 11. PRINT TOP 15 UP TERMS
top_up = up_sig.sort_values(
    "Adjusted P-value"
).head(15)

print("\n================================")
print("TOP 15 UPREGULATED GO TERMS")
print("================================\n")

print(
    top_up[
        [
            "Term",
            "Adjusted P-value",
            "Odds Ratio",
            "Combined Score"
        ]
    ].to_string(index=False)
)

# 12. PRINT TOP 15 DOWN TERMS
top_down = down_sig.sort_values(
    "Adjusted P-value"
).head(15)

print("\n================================")
print("TOP 15 DOWNREGULATED GO TERMS")
print("================================\n")

print(
    top_down[
        [
            "Term",
            "Adjusted P-value",
            "Odds Ratio",
            "Combined Score"
        ]
    ].to_string(index=False)
)

# 13. FUNCTION TO DRAW GO BARPLOT
def plot_go_terms(df, title, output_file):

    if len(df) == 0:
        print(f"No significant terms for {title}")
        return

    plot_df = (
        df.sort_values(
            "Adjusted P-value"
        )
        .head(15)
        .copy()
    )

    plot_df["minus_log10_padj"] = -np.log10(
        plot_df["Adjusted P-value"].clip(
            lower=1e-300
        )
    )
    plot_df = plot_df.iloc[::-1]

    plt.figure(
        figsize=(10, 8)
    )

    plt.barh(
        plot_df["Term"],
        plot_df["minus_log10_padj"]
    )

    plt.xlabel(
        "-log10 adjusted p-value"
    )

    plt.ylabel(
        "GO Biological Process"
    )

    plt.title(
        title
    )

    plt.tight_layout()

    plt.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

# 14. SAVE GO FIGURES
plot_go_terms(
    up_sig,
    "GO Enrichment of Upregulated Genes",
    "results/figures/go_upregulated.png"
)

plot_go_terms(
    down_sig,
    "GO Enrichment of Downregulated Genes",
    "results/figures/go_downregulated.png"
)


print("\nGO enrichment analysis completed.")