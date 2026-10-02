import pandas as pd
from pathlib import Path

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# 1. PATHS
COUNTS_FILE = Path(
    "data/processed/filtered_counts.csv"
)

METADATA_FILE = Path(
    "data/processed/metadata.csv"
)

TABLE_DIR = Path(
    "results/tables"
)

TABLE_DIR.mkdir(
    parents=True,
    exist_ok=True
)

# 2. ANALYSIS PARAMETERS
ALPHA = 0.05
LOG2FC_THRESHOLD = 1.0

# 3. LOAD DATA
counts = pd.read_csv(
    COUNTS_FILE,
    index_col=0
)

metadata = pd.read_csv(
    METADATA_FILE,
    index_col=0
)

counts = counts.astype(int)


print(
    "Count matrix:"
)

print(
    f"{counts.shape[0]} samples x "
    f"{counts.shape[1]} genes"
)

# 4. VALIDATE DATA
if list(counts.index) != list(metadata.index):
    raise ValueError(
        "Counts and metadata sample order do not match."
    )

if (counts < 0).any().any():
    raise ValueError(
        "Negative counts detected."
    )

print(
    "Input validation: PASSED"
)

# 5. CREATE DESEQ DATASET
dds = DeseqDataSet(
    counts=counts,
    metadata=metadata,
    design="~condition",
    refit_cooks=True,
    n_cpus=1
)


# 6. FIT MODEL
print(
    "\nRunning PyDESeq2..."
)

dds.deseq2()

print(
    "Model fitting completed."
)

# 7. DIFFERENTIAL EXPRESSION TEST
stats = DeseqStats(
    dds,
    contrast=[
        "condition",
        "mRNA_LNP",
        "control"
    ],
    alpha=ALPHA,
    n_cpus=1
)

stats.summary()

results = (
    stats.results_df.copy()
)

# 8. SAVE ALL RESULTS
results.to_csv(
    TABLE_DIR /
    "pydeseq2_all_results.csv"
)

# 9. DEFINE SIGNIFICANT DEGs
significant = results[
    (results["padj"] < ALPHA)
    &
    (
        results[
            "log2FoldChange"
        ].abs()
        >= LOG2FC_THRESHOLD
    )
].copy()


upregulated = significant[
    significant[
        "log2FoldChange"
    ] >= LOG2FC_THRESHOLD
].copy()


downregulated = significant[
    significant[
        "log2FoldChange"
    ] <= -LOG2FC_THRESHOLD
].copy()

# 10. SAVE DEG TABLES
significant.to_csv(
    TABLE_DIR /
    "pydeseq2_significant_DEGs.csv"
)

upregulated.to_csv(
    TABLE_DIR /
    "pydeseq2_upregulated_DEGs.csv"
)

downregulated.to_csv(
    TABLE_DIR /
    "pydeseq2_downregulated_DEGs.csv"
)

# 11. eGFP SANITY CHECK
print(
    "\neGFP differential expression:"
)

if "eGFP" in results.index:

    egfp_result = results.loc[
        "eGFP",
        [
            "baseMean",
            "log2FoldChange",
            "pvalue",
            "padj"
        ]
    ]

    print(
        egfp_result
    )

else:

    print(
        "eGFP not found."
    )

# 12. TOP DEGs
top_up = (
    upregulated
    .sort_values(
        "log2FoldChange",
        ascending=False
    )
    .head(10)
)

top_down = (
    downregulated
    .sort_values(
        "log2FoldChange",
        ascending=True
    )
    .head(10)
)


print(
    "\nTop 10 upregulated genes:"
)

print(
    top_up[
        [
            "baseMean",
            "log2FoldChange",
            "padj"
        ]
    ]
)


print(
    "\nTop 10 downregulated genes:"
)

print(
    top_down[
        [
            "baseMean",
            "log2FoldChange",
            "padj"
        ]
    ]
)

# 13. FINAL SUMMARY
print(
    "\n================================"
)

print(
    "PYDESEQ2 DEG SUMMARY"
)

print(
    "================================"
)

print(
    f"Genes tested: "
    f"{len(results)}"
)

print(
    f"Significant DEGs: "
    f"{len(significant)}"
)

print(
    f"Upregulated: "
    f"{len(upregulated)}"
)

print(
    f"Downregulated: "
    f"{len(downregulated)}"
)

print(
    "\nThresholds:"
)

print(
    f"Adjusted p-value < {ALPHA}"
)

print(
    f"|log2FoldChange| >= "
    f"{LOG2FC_THRESHOLD}"
)

print(
    "\nPyDESeq2 analysis completed."
)