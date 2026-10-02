from pathlib import Path

import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt

# 1. PATHS
PYDESEQ2_FILE = Path(
    "results/tables/pydeseq2_all_results.csv"
)

EDGER_FILE = Path(
    "results/tables/edger_all_results.csv"
)

ANNOTATION_FILE = Path(
    "results/tables/pydeseq2_annotated_results.csv"
)

TABLE_DIR = Path(
    "results/tables"
)

FIGURE_DIR = Path(
    "results/figures"
)

TABLE_DIR.mkdir(
    parents=True,
    exist_ok=True
)

FIGURE_DIR.mkdir(
    parents=True,
    exist_ok=True
)

# 2. ANALYSIS PARAMETERS
SIGNIFICANCE_THRESHOLD = 0.05
LOG2FC_THRESHOLD = 1.0

# 3. LOAD RESULTS
print(
    "Loading PyDESeq2 and edgeR results..."
)

py = pd.read_csv(
    PYDESEQ2_FILE,
    index_col=0
)

edge = pd.read_csv(
    EDGER_FILE,
    index_col=0
)

annot = pd.read_csv(
    ANNOTATION_FILE
)


print(
    f"PyDESeq2 genes: {len(py)}"
)

print(
    f"edgeR genes: {len(edge)}"
)

# 4. CHECK GENE SET CONSISTENCY
py_gene_set = set(
    py.index
)

edge_gene_set = set(
    edge.index
)

if py_gene_set != edge_gene_set:

    missing_in_edge = (
        py_gene_set
        - edge_gene_set
    )

    missing_in_py = (
        edge_gene_set
        - py_gene_set
    )

    raise ValueError(
        "PyDESeq2 and edgeR did not test "
        "the same gene set.\n"
        f"Missing in edgeR: "
        f"{len(missing_in_edge)}\n"
        f"Missing in PyDESeq2: "
        f"{len(missing_in_py)}"
    )


print(
    "Gene set consistency check: PASSED"
)

# 5. CREATE GENE SYMBOL MAPPING
gene_symbol_map = (
    annot[
        [
            "original_gene_id",
            "gene_symbol"
        ]
    ]
    .drop_duplicates(
        "original_gene_id"
    )
    .set_index(
        "original_gene_id"
    )["gene_symbol"]
    .to_dict()
)

# 6. DEFINE SIGNIFICANT DEGs
py_sig = (
    (py["padj"] < SIGNIFICANCE_THRESHOLD)
    &
    (
        py["log2FoldChange"].abs()
        >= LOG2FC_THRESHOLD
    )
)

edge_sig = (
    (edge["FDR"] < SIGNIFICANCE_THRESHOLD)
    &
    (
        edge["logFC"].abs()
        >= LOG2FC_THRESHOLD
    )
)


py_genes = set(
    py.index[
        py_sig
    ]
)

edge_genes = set(
    edge.index[
        edge_sig
    ]
)

# 7. IDENTIFY SHARED AND METHOD-SPECIFIC DEGs
shared_genes = (
    py_genes
    & edge_genes
)

py_only = sorted(
    py_genes
    - edge_genes
)

edge_only = sorted(
    edge_genes
    - py_genes
)

discordant_genes = (
    py_only
    + edge_only
)

# 8. BUILD DISCORDANT GENE TABLE
comparison = pd.DataFrame(
    index=discordant_genes
)

comparison.index.name = (
    "gene_id"
)


comparison["gene_symbol"] = [
    gene_symbol_map.get(
        gene
    )
    for gene in comparison.index
]


comparison[
    "PyDESeq2_log2FC"
] = py.loc[
    comparison.index,
    "log2FoldChange"
]


comparison[
    "PyDESeq2_padj"
] = py.loc[
    comparison.index,
    "padj"
]


comparison[
    "edgeR_logFC"
] = edge.loc[
    comparison.index,
    "logFC"
]


comparison[
    "edgeR_FDR"
] = edge.loc[
    comparison.index,
    "FDR"
]


comparison[
    "classification"
] = [
    (
        "PyDESeq2-only"
        if gene in py_only
        else "edgeR-only"
    )
    for gene in comparison.index
]

# 9. DETERMINE REASON FOR METHOD DISAGREEMENT
def determine_reason(row):

    if (
        row["classification"]
        == "PyDESeq2-only"
    ):

        fail_significance = (
            pd.isna(
                row["edgeR_FDR"]
            )
            or
            row["edgeR_FDR"]
            >= SIGNIFICANCE_THRESHOLD
        )

        fail_fold_change = (
            pd.isna(
                row["edgeR_logFC"]
            )
            or
            abs(
                row["edgeR_logFC"]
            )
            < LOG2FC_THRESHOLD
        )

    else:

        fail_significance = (
            pd.isna(
                row["PyDESeq2_padj"]
            )
            or
            row["PyDESeq2_padj"]
            >= SIGNIFICANCE_THRESHOLD
        )

        fail_fold_change = (
            pd.isna(
                row["PyDESeq2_log2FC"]
            )
            or
            abs(
                row["PyDESeq2_log2FC"]
            )
            < LOG2FC_THRESHOLD
        )


    if (
        fail_significance
        and fail_fold_change
    ):

        return (
            "Fails both significance "
            "and FC thresholds"
        )

    elif fail_significance:

        return (
            "Fails significance threshold"
        )

    elif fail_fold_change:

        return (
            "Fails fold-change threshold"
        )

    else:

        return (
            "Other / numerical boundary case"
        )


comparison[
    "reason"
] = comparison.apply(
    determine_reason,
    axis=1
)

# 10. CALCULATE DISTANCE FROM FC THRESHOLD
comparison[
    "PyDESeq2_absFC"
] = comparison[
    "PyDESeq2_log2FC"
].abs()


comparison[
    "edgeR_absFC"
] = comparison[
    "edgeR_logFC"
].abs()


comparison[
    "PyDESeq2_distance_from_FC_threshold"
] = (
    comparison[
        "PyDESeq2_absFC"
    ]
    - LOG2FC_THRESHOLD
)


comparison[
    "edgeR_distance_from_FC_threshold"
] = (
    comparison[
        "edgeR_absFC"
    ]
    - LOG2FC_THRESHOLD
)

# 11. PRINT SUMMARY
print(
    "\n================================"
)

print(
    "DISCORDANT DEG ANALYSIS"
)

print(
    "================================"
)

print(
    f"Shared DEGs: "
    f"{len(shared_genes)}"
)

print(
    f"PyDESeq2-only DEGs: "
    f"{len(py_only)}"
)

print(
    f"edgeR-only DEGs: "
    f"{len(edge_only)}"
)

print(
    f"Total discordant genes: "
    f"{len(comparison)}"
)

# 12. SUMMARIZE REASONS
reason_counts = (
    comparison[
        "reason"
    ]
    .value_counts()
)


print(
    "\nReasons for disagreement:"
)

print(
    reason_counts
)

# 13. PRINT DISCORDANT TABLE
display_columns = [
    "gene_symbol",
    "classification",
    "PyDESeq2_log2FC",
    "PyDESeq2_padj",
    "edgeR_logFC",
    "edgeR_FDR",
    "reason"
]


print(
    "\nDiscordant genes:\n"
)

print(
    comparison[
        display_columns
    ]
    .sort_values(
        [
            "classification",
            "gene_symbol"
        ],
        na_position="last"
    )
    .to_string()
)

# 14. SAVE DISCORDANT GENE TABLE
DISCORDANT_FILE = (
    TABLE_DIR
    / "discordant_DEGs_pydeseq2_vs_edger.csv"
)

comparison.to_csv(
    DISCORDANT_FILE
)

# 15. SAVE DISAGREEMENT SUMMARY TABLE
reason_summary = (
    reason_counts
    .rename_axis(
        "reason"
    )
    .reset_index(
        name="gene_count"
    )
)

reason_summary.to_csv(
    TABLE_DIR
    / "discordant_DEG_reason_summary.csv",
    index=False
)

# 16. PLOT REASONS FOR DISAGREEMENT
plt.figure(
    figsize=(9, 6)
)

plt.bar(
    reason_counts.index,
    reason_counts.values
)

plt.ylabel(
    "Number of genes"
)

plt.xlabel(
    "Reason for method disagreement"
)

plt.title(
    "Why PyDESeq2 and edgeR "
    "Identified Different DEGs"
)

plt.xticks(
    rotation=20,
    ha="right"
)


for i, value in enumerate(
    reason_counts.values
):

    plt.text(
        i,
        value + 0.3,
        str(value),
        ha="center"
    )


plt.tight_layout()


DISCORDANT_FIGURE = (
    FIGURE_DIR
    / "discordant_deg_reasons.png"
)

plt.savefig(
    DISCORDANT_FIGURE,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

# 17. FINAL SUMMARY
print(
    "\n================================"
)

print(
    "DISCORDANT ANALYSIS COMPLETED"
)

print(
    "================================"
)

print(
    f"Discordant genes table:\n"
    f"{DISCORDANT_FILE}"
)

print(
    "\nReason summary table:\n"
    f"{TABLE_DIR / 'discordant_DEG_reason_summary.csv'}"
)

print(
    "\nFigure:\n"
    f"{DISCORDANT_FIGURE}"
)