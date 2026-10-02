import pandas as pd
import numpy as np

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt

# 1. LOAD FILTERED COUNTS
counts = pd.read_csv(
    "data/processed/filtered_counts.csv",
    index_col=0
)

# counts currently = samples x genes

# 2. LOAD ANNOTATION TABLE
annotated = pd.read_csv(
    "results/tables/pydeseq2_annotated_results.csv"
)

# Create mapping:
# Ensembl/original ID -> gene symbol
gene_map = annotated.set_index(
    "original_gene_id"
)["gene_symbol"].to_dict()

# 3. SELECT BIOLOGICALLY RELEVANT GENES
selected_symbols = [
    "Isg15",
    "Ifi44",
    "Ifit3b",
    "Ifi27l2b",
    "Tgtp1",
    "Cxcl9",
    "Ccl7",
    "Saa3",
    "Ly6a",
    "Cdkn1a"
]

# 4. FIND MATCHING GENE IDS
selected_gene_ids = []

for gene_id, symbol in gene_map.items():

    if symbol in selected_symbols:
        selected_gene_ids.append(gene_id)


print("Selected genes found:")
print(len(selected_gene_ids))

# 5. EXTRACT COUNTS
selected_counts = counts[
    selected_gene_ids
].copy()

# 6. RENAME COLUMNS TO GENE SYMBOLS
selected_counts = selected_counts.rename(
    columns=gene_map
)

# 7. LOG TRANSFORM
log_counts = np.log2(
    selected_counts + 1
)

# 8. Z-SCORE EACH GENE
z_scores = (
    log_counts - log_counts.mean(axis=0)
) / log_counts.std(axis=0)

# 9. TRANSPOSE FOR HEATMAP
heatmap_data = z_scores.T

# 10. DRAW HEATMAP
plt.figure(figsize=(10, 7))

plt.imshow(
    heatmap_data,
    aspect="auto",
    interpolation="nearest"
)

plt.colorbar(
    label="Z-score"
)

plt.xticks(
    range(len(heatmap_data.columns)),
    heatmap_data.columns,
    rotation=45,
    ha="right"
)

plt.yticks(
    range(len(heatmap_data.index)),
    heatmap_data.index
)

plt.title(
    "Expression of Selected Immune-Response Genes"
)

plt.xlabel(
    "Samples"
)

plt.ylabel(
    "Genes"
)

plt.tight_layout()

# 11. SAVE
plt.savefig(
    "results/figures/selected_genes_heatmap.png",
    dpi=300
)

plt.close()

print("\nHeatmap saved:")
print(
    "results/figures/selected_genes_heatmap.png"
)