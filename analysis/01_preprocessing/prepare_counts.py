import pandas as pd
from pathlib import Path


# 1. PROJECT PATHS
RAW_FILE = Path(
    "data/raw/GSE331154_mmuegfp_liver_counts.txt"
)

PROCESSED_DIR = Path(
    "data/processed"
)

FILTERED_COUNTS_FILE = (
    PROCESSED_DIR / "filtered_counts.csv"
)

METADATA_FILE = (
    PROCESSED_DIR / "metadata.csv"
)

# 2. EXPERIMENTAL DESIGN
CONTROL_SAMPLES = [
    "saline1",
    "saline2",
    "saline3",
    "saline4",
    "saline5",
]

TREATMENT_SAMPLES = [
    "mRNA-eGFP1",
    "mRNA-eGFP2",
    "mRNA-eGFP3",
    "mRNA-eGFP4",
    "mRNA-eGFP5",
]

SAMPLE_COLUMNS = (
    CONTROL_SAMPLES
    + TREATMENT_SAMPLES
)


# 3. FILTERING PARAMETERS
MIN_COUNT = 10
MIN_SAMPLES = 3

# 4. CHECK INPUT FILE

if not RAW_FILE.exists():
    raise FileNotFoundError(
        f"Raw count file not found: {RAW_FILE}"
    )

# 5. READ featureCounts OUTPUT
print(
    "Reading raw featureCounts table..."
)

raw = pd.read_csv(
    RAW_FILE,
    sep="\t"
)

print(
    f"Raw dataset shape: {raw.shape}"
)

# 6. VALIDATE REQUIRED COLUMNS
required_columns = (
    ["Geneid"]
    + SAMPLE_COLUMNS
)

missing_columns = [
    column
    for column in required_columns
    if column not in raw.columns
]

if missing_columns:
    raise ValueError(
        "Missing required columns: "
        + ", ".join(missing_columns)
    )

print(
    "Required sample columns: PASSED"
)

# 7. CHECK DUPLICATED GENE IDs
duplicate_gene_ids = (
    raw["Geneid"].duplicated().sum()
)

print(
    f"Duplicated Gene IDs: "
    f"{duplicate_gene_ids}"
)

if duplicate_gene_ids > 0:
    raise ValueError(
        "Duplicated Gene IDs detected. "
        "Check the input file before continuing."
    )


# 8. CREATE COUNT MATRIX
counts = raw[
    ["Geneid"] + SAMPLE_COLUMNS
].copy()

counts = counts.set_index(
    "Geneid"
)

# PyDESeq2 expects:
# samples x genes

counts = counts.T


print(
    "\nCount matrix before filtering:"
)

print(
    f"{counts.shape[0]} samples x "
    f"{counts.shape[1]} genes"
)

# 9. CHECK COUNTS
if (counts < 0).any().any():
    raise ValueError(
        "Negative count values detected."
    )

counts = counts.astype(int)

# 10. FILTER LOW-COUNT GENES
keep = (
    (counts >= MIN_COUNT)
    .sum(axis=0)
    >= MIN_SAMPLES
)

filtered_counts = counts.loc[
    :,
    keep
].copy()


genes_before = counts.shape[1]

genes_after = (
    filtered_counts.shape[1]
)

genes_removed = (
    genes_before - genes_after
)


print(
    "\nFiltering rule:"
)

print(
    f"Keep genes with count >= "
    f"{MIN_COUNT} in at least "
    f"{MIN_SAMPLES} samples"
)

print(
    f"\nGenes before filtering: "
    f"{genes_before}"
)

print(
    f"Genes after filtering: "
    f"{genes_after}"
)

print(
    f"Genes removed: "
    f"{genes_removed}"
)

# 11. CREATE METADATA
metadata = pd.DataFrame(
    {
        "condition": (
            ["control"] * len(
                CONTROL_SAMPLES
            )
            +
            ["mRNA_LNP"] * len(
                TREATMENT_SAMPLES
            )
        )
    },
    index=SAMPLE_COLUMNS
)

metadata.index.name = "sample"

# 12. VERIFY SAMPLE ORDER
if list(filtered_counts.index) != SAMPLE_COLUMNS:
    raise ValueError(
        "Sample order in count matrix "
        "does not match metadata."
    )

print(
    "\nSample order check: PASSED"
)

# 13. CREATE OUTPUT DIRECTORY
PROCESSED_DIR.mkdir(
    parents=True,
    exist_ok=True
)

# 14. SAVE PROCESSED FILES
filtered_counts.to_csv(
    FILTERED_COUNTS_FILE
)

metadata.to_csv(
    METADATA_FILE
)

# 15. FINAL SUMMARY
print(
    "\n================================"
)

print(
    "PREPROCESSING COMPLETED"
)

print(
    "================================"
)

print(
    f"Filtered count matrix:"
)

print(
    FILTERED_COUNTS_FILE
)

print(
    f"\nMetadata:"
)

print(
    METADATA_FILE
)

print(
    f"\nFinal matrix shape: "
    f"{filtered_counts.shape[0]} "
    f"samples x "
    f"{filtered_counts.shape[1]} genes"
)

print(
    "\nMetadata:"
)

print(
    metadata
)