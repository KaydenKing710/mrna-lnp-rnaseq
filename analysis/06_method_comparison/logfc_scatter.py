import pandas as pd
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

# 1. LOAD RESULTS
pydeseq2 = pd.read_csv(
    "results/tables/pydeseq2_all_results.csv",
    index_col=0
)

edger = pd.read_csv(
    "results/tables/edger_all_results.csv",
    index_col=0
)

# 2. MATCH GENES
comparison = pd.DataFrame({
    "PyDESeq2": pydeseq2["log2FoldChange"],
    "edgeR": edger["logFC"]
}).dropna()

# 3. CALCULATE CORRELATION
pearson = comparison["PyDESeq2"].corr(
    comparison["edgeR"],
    method="pearson"
)

# 4. CREATE SCATTER PLOT
plt.figure(figsize=(8, 8))

plt.scatter(
    comparison["PyDESeq2"],
    comparison["edgeR"],
    s=12,
    alpha=0.5
)

# 5. ADD PERFECT AGREEMENT LINE
minimum = min(
    comparison["PyDESeq2"].min(),
    comparison["edgeR"].min()
)

maximum = max(
    comparison["PyDESeq2"].max(),
    comparison["edgeR"].max()
)

plt.plot(
    [minimum, maximum],
    [minimum, maximum],
    linestyle="--"
)

# 6. LABEL eGFP
if "eGFP" in comparison.index:

    egfp = comparison.loc["eGFP"]

    plt.annotate(
        "eGFP",
        (egfp["PyDESeq2"], egfp["edgeR"]),
        xytext=(5, 5),
        textcoords="offset points"
    )


# 7. LABEL FIGURE
plt.xlabel(
    "PyDESeq2 log2 Fold Change"
)

plt.ylabel(
    "edgeR log2 Fold Change"
)

plt.title(
    f"Fold-Change Agreement: PyDESeq2 vs edgeR\n"
    f"Pearson r = {pearson:.4f}"
)

plt.tight_layout()
# 8. SAVE
plt.savefig(
    "results/figures/pydeseq2_vs_edger_logfc.png",
    dpi=300
)

plt.close()

print(
    "Comparison plot saved to:"
)

print(
    "results/figures/pydeseq2_vs_edger_logfc.png"
)