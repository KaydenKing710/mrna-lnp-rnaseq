import pandas as pd
import numpy as np
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# 1. PATHS
COUNTS_FILE = Path(
    "data/processed/filtered_counts.csv"
)

METADATA_FILE = Path(
    "data/processed/metadata.csv"
)

FIGURE_DIR = Path(
    "results/figures"
)

TABLE_DIR = Path(
    "results/tables"
)

FIGURE_DIR.mkdir(
    parents=True,
    exist_ok=True
)

TABLE_DIR.mkdir(
    parents=True,
    exist_ok=True
)
# 2. LOAD DATA
counts = pd.read_csv(
    COUNTS_FILE,
    index_col=0
)

metadata = pd.read_csv(
    METADATA_FILE,
    index_col=0
)

print(
    "Count matrix:"
)

print(
    f"{counts.shape[0]} samples x "
    f"{counts.shape[1]} genes"
)

print(
    "\nMetadata:"
)

print(metadata)

# 3. VALIDATE SAMPLE MATCHING

if list(counts.index) != list(metadata.index):
    raise ValueError(
        "Sample names/order in counts and metadata "
        "do not match."
    )

print(
    "\nSample matching check: PASSED"
)

# 4. LIBRARY SIZE QC

library_sizes = counts.sum(
    axis=1
)

library_sizes.to_csv(
    TABLE_DIR / "library_sizes.csv",
    header=["library_size"]
)

print(
    "\nLibrary sizes:"
)

print(
    library_sizes
)

# 5. CPM NORMALIZATION FOR PCA
cpm = counts.div(
    library_sizes,
    axis=0
) * 1_000_000

log_cpm = np.log2(
    cpm + 1
)

# 6. PCA
pca = PCA(
    n_components=2
)

pca_result = pca.fit_transform(
    log_cpm
)

pca_df = pd.DataFrame(
    pca_result,
    columns=["PC1", "PC2"],
    index=counts.index
)

pca_df["condition"] = (
    metadata["condition"]
)

pca_df.to_csv(
    TABLE_DIR / "pca_coordinates.csv"
)


pc1_variance = (
    pca.explained_variance_ratio_[0]
)

pc2_variance = (
    pca.explained_variance_ratio_[1]
)


print(
    "\nPCA explained variance:"
)

print(
    f"PC1: {pc1_variance * 100:.2f}%"
)

print(
    f"PC2: {pc2_variance * 100:.2f}%"
)

# 7. PCA FIGURE
plt.figure(
    figsize=(8, 6)
)

for condition in pca_df[
    "condition"
].unique():

    subset = pca_df[
        pca_df["condition"] == condition
    ]

    plt.scatter(
        subset["PC1"],
        subset["PC2"],
        label=condition,
        s=80
    )

    for sample in subset.index:

        plt.annotate(
            sample,
            (
                subset.loc[sample, "PC1"],
                subset.loc[sample, "PC2"]
            ),
            xytext=(3, 3),
            textcoords="offset points",
            fontsize=8
        )


plt.xlabel(
    f"PC1 ({pc1_variance * 100:.1f}%)"
)

plt.ylabel(
    f"PC2 ({pc2_variance * 100:.1f}%)"
)

plt.title(
    "PCA of mRNA-LNP vs Saline RNA-seq Samples"
)

plt.legend()

plt.tight_layout()

plt.savefig(
    FIGURE_DIR / "pca_plot.png",
    dpi=300
)

plt.close()

# 8. eGFP SANITY CHECK

if "eGFP" not in counts.columns:
    raise ValueError(
        "eGFP was not found in the filtered count matrix."
    )

egfp_counts = counts[
    "eGFP"
].copy()

egfp_counts.to_csv(
    TABLE_DIR / "egfp_counts.csv",
    header=["eGFP_raw_count"]
)

print(
    "\neGFP raw counts:"
)

print(
    egfp_counts
)

# 9. eGFP GROUP SUMMARY
egfp_df = pd.DataFrame({
    "eGFP_count": egfp_counts,
    "condition": metadata["condition"]
})

egfp_summary = (
    egfp_df.groupby("condition")[
        "eGFP_count"
    ]
    .agg(
        ["mean", "median", "min", "max"]
    )
)

egfp_summary.to_csv(
    TABLE_DIR / "egfp_group_summary.csv"
)

print(
    "\neGFP group summary:"
)

print(
    egfp_summary
)

# 10. FINAL SUMMARY
print(
    "\n================================"
)

print(
    "QUALITY CONTROL COMPLETED"
)

print(
    "================================"
)

print(
    f"PCA figure: "
    f"{FIGURE_DIR / 'pca_plot.png'}"
)

print(
    f"PC1 variance: "
    f"{pc1_variance * 100:.2f}%"
)

print(
    f"PC2 variance: "
    f"{pc2_variance * 100:.2f}%"
)

print(
    "\nQC tables saved in:"
)

print(
    TABLE_DIR
)