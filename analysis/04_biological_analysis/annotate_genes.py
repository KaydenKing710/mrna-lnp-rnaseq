import pandas as pd
import mygene
from pathlib import Path

# 1. PATHS
RESULTS_FILE = Path(
    "results/tables/pydeseq2_all_results.csv"
)

TABLE_DIR = Path(
    "results/tables"
)

TABLE_DIR.mkdir(
    parents=True,
    exist_ok=True
)
# 2. PARAMETERS
ALPHA = 0.05
LOG2FC_THRESHOLD = 1.0

# 3. LOAD PYDESEQ2 RESULTS
results = pd.read_csv(
    RESULTS_FILE,
    index_col=0
)

print(
    f"Genes loaded: {len(results)}"
)

# 4. PREPARE ENSEMBL IDs
results["original_gene_id"] = (
    results.index
)

results["ensembl_id"] = (
    results["original_gene_id"]
    .str.split(".")
    .str[0]
)


ensembl_mask = (
    results["ensembl_id"]
    .str.startswith("ENSMUSG")
)

ensembl_ids = (
    results.loc[
        ensembl_mask,
        "ensembl_id"
    ]
    .drop_duplicates()
    .tolist()
)

print(
    f"Mouse Ensembl IDs to annotate: "
    f"{len(ensembl_ids)}"
)

# 5. QUERY MYGENE
mg = mygene.MyGeneInfo()

annotations = mg.querymany(
    ensembl_ids,
    scopes="ensembl.gene",
    fields="symbol,name",
    species="mouse",
    as_dataframe=True
)
# 6. CLEAN ANNOTATIONS
annotations = (
    annotations
    .reset_index()
    .rename(
        columns={
            "query": "ensembl_id",
            "symbol": "gene_symbol",
            "name": "gene_name"
        }
    )
)

annotations = (
    annotations[
        [
            "ensembl_id",
            "gene_symbol",
            "gene_name"
        ]
    ]
    .drop_duplicates(
        subset="ensembl_id"
    )
)

# 7. MERGE WITH DEG RESULTS
annotated = (
    results
    .reset_index(drop=True)
    .merge(
        annotations,
        on="ensembl_id",
        how="left"
    )
)

# 8. PRESERVE eGFP
egfp_mask = (
    annotated["original_gene_id"]
    == "eGFP"
)

annotated.loc[
    egfp_mask,
    "gene_symbol"
] = "eGFP"

annotated.loc[
    egfp_mask,
    "gene_name"
] = (
    "enhanced green fluorescent protein"
)
# 9. SAVE FULL ANNOTATED RESULTS
ANNOTATED_FILE = (
    TABLE_DIR /
    "pydeseq2_annotated_results.csv"
)

annotated.to_csv(
    ANNOTATED_FILE,
    index=False
)


mapped_count = (
    annotated[
        "gene_symbol"
    ]
    .notna()
    .sum()
)

print(
    f"\nGenes mapped to symbols: "
    f"{mapped_count}"
)

# 10. SELECT SIGNIFICANT DEGs
significant = annotated[
    (annotated["padj"] < ALPHA)
    &
    (
        annotated[
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

# 11. TOP ANNOTATED GENES
top_up = (
    upregulated
    .dropna(
        subset=["gene_symbol"]
    )
    .sort_values(
        "log2FoldChange",
        ascending=False
    )
    .head(20)
)


top_down = (
    downregulated
    .dropna(
        subset=["gene_symbol"]
    )
    .sort_values(
        "log2FoldChange",
        ascending=True
    )
    .head(20)
)

# 12. SAVE TOP GENE TABLES
top_up.to_csv(
    TABLE_DIR /
    "top20_upregulated.csv",
    index=False
)

top_down.to_csv(
    TABLE_DIR /
    "top20_downregulated.csv",
    index=False
)

# 13. PRINT TOP RESULTS
display_columns = [
    "gene_symbol",
    "gene_name",
    "baseMean",
    "log2FoldChange",
    "padj"
]


print(
    "\n================================"
)

print(
    "TOP 20 UPREGULATED GENES"
)

print(
    "================================\n"
)

print(
    top_up[
        display_columns
    ].to_string(
        index=False
    )
)


print(
    "\n================================"
)

print(
    "TOP 20 DOWNREGULATED GENES"
)

print(
    "================================\n"
)

print(
    top_down[
        display_columns
    ].to_string(
        index=False
    )
)

# 14. FINAL SUMMARY
print(
    "\n================================"
)

print(
    "GENE ANNOTATION COMPLETED"
)

print(
    "================================"
)

print(
    f"Annotated results: "
    f"{ANNOTATED_FILE}"
)

print(
    f"Significant annotated DEGs: "
    f"{len(significant)}"
)