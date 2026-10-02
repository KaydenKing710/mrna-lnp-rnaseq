import pandas as pd

# 1. LOAD RESULTS
pydeseq2 = pd.read_csv(
    "results/tables/pydeseq2_all_results.csv",
    index_col=0
)

edger = pd.read_csv(
    "results/tables/edger_all_results.csv",
    index_col=0
)


print("PyDESeq2 genes:")
print(len(pydeseq2))

print("\nedgeR genes:")
print(len(edger))


# 2. DEFINE SIGNIFICANT DEGs
# SAME THRESHOLDS FOR BOTH METHODS
py_sig = pydeseq2[
    (pydeseq2["padj"] < 0.05) &
    (pydeseq2["log2FoldChange"].abs() >= 1)
].copy()


edge_sig = edger[
    (edger["FDR"] < 0.05) &
    (edger["logFC"].abs() >= 1)
].copy()


print("\nPyDESeq2 significant DEGs:")
print(len(py_sig))

print("\nedgeR significant DEGs:")
print(len(edge_sig))

# 3. CREATE DEG SETS
py_genes = set(py_sig.index)

edge_genes = set(edge_sig.index)


shared = py_genes & edge_genes

py_only = py_genes - edge_genes

edge_only = edge_genes - py_genes

union = py_genes | edge_genes

# 4. JACCARD SIMILARITY
jaccard = len(shared) / len(union)

# 5. PRINT OVERLAP SUMMARY
print("\n================================")
print("DEG SET COMPARISON")
print("================================")

print("Shared DEGs:")
print(len(shared))

print("\nPyDESeq2-only DEGs:")
print(len(py_only))

print("\nedgeR-only DEGs:")
print(len(edge_only))

print("\nUnion of DEGs:")
print(len(union))

print("\nJaccard similarity:")
print(round(jaccard, 4))

# 6. OVERLAP PERCENTAGES
py_overlap = len(shared) / len(py_genes) * 100

edge_overlap = len(shared) / len(edge_genes) * 100

print("\nPercent of PyDESeq2 DEGs also found by edgeR:")
print(round(py_overlap, 2), "%")

print("\nPercent of edgeR DEGs also found by PyDESeq2:")
print(round(edge_overlap, 2), "%")

# 7. SAVE DEG LISTS
pd.DataFrame(
    sorted(shared),
    columns=["gene_id"]
).to_csv(
    "results/tables/shared_DEGs.csv",
    index=False
)

pd.DataFrame(
    sorted(py_only),
    columns=["gene_id"]
).to_csv(
    "results/tables/pydeseq2_only_DEGs.csv",
    index=False
)

pd.DataFrame(
    sorted(edge_only),
    columns=["gene_id"]
).to_csv(
    "results/tables/edger_only_DEGs.csv",
    index=False
)

# 8. COMPARE LOG2 FOLD CHANGES
# FOR ALL GENES
comparison = pd.DataFrame({
    "pydeseq2_log2FC": pydeseq2["log2FoldChange"],
    "edger_log2FC": edger["logFC"]
}).dropna()


pearson = comparison[
    "pydeseq2_log2FC"
].corr(
    comparison["edger_log2FC"],
    method="pearson"
)

spearman = comparison[
    "pydeseq2_log2FC"
].corr(
    comparison["edger_log2FC"],
    method="spearman"
)

print("\n================================")
print("LOG2 FOLD-CHANGE AGREEMENT")
print("================================")

print("Pearson correlation:")
print(round(pearson, 4))

print("\nSpearman correlation:")
print(round(spearman, 4))

# 9. DIRECTION AGREEMENT AMONG SHARED DEGs
shared_list = list(shared)

shared_fc = pd.DataFrame({
    "pydeseq2_log2FC":
        pydeseq2.loc[shared_list, "log2FoldChange"],

    "edger_log2FC":
        edger.loc[shared_list, "logFC"]
})


same_direction = (
    shared_fc["pydeseq2_log2FC"] *
    shared_fc["edger_log2FC"] > 0
)

direction_agreement = (
    same_direction.mean() * 100
)

print("\nDirection agreement among shared DEGs:")
print(round(direction_agreement, 2), "%")

# 10. SAVE FULL COMPARISON TABLE
comparison.to_csv(
    "results/tables/log2fc_comparison.csv"
)

print("\nComparison files saved successfully.")