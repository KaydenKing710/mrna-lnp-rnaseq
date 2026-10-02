import pandas as pd
import numpy as np

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt

# 1. LOAD PYDESEQ2 RESULTS
results = pd.read_csv(
    "results/tables/pydeseq2_all_results.csv",
    index_col=0
)

print("Genes loaded:")
print(len(results))

# 2. REMOVE GENES WITHOUT PADJ
plot_df = results.dropna(
    subset=["padj", "log2FoldChange"]
).copy()

# 3. DEFINE DEG STATUS
plot_df["status"] = "Not significant"

plot_df.loc[
    (plot_df["padj"] < 0.05) &
    (plot_df["log2FoldChange"] >= 1),
    "status"
] = "Upregulated"

plot_df.loc[
    (plot_df["padj"] < 0.05) &
    (plot_df["log2FoldChange"] <= -1),
    "status"
] = "Downregulated"

# 4. CALCULATE -LOG10(PADJ)
# Some extremely small adjusted p-values may be stored as zero due to numerical precision.
plot_df["padj_safe"] = plot_df["padj"].clip(
    lower=1e-50
)

plot_df["minus_log10_padj"] = -np.log10(
    plot_df["padj_safe"]
)

# 5. CREATE VOLCANO PLOT
plt.figure(figsize=(9, 7))

for status in [
    "Not significant",
    "Upregulated",
    "Downregulated"
]:

    subset = plot_df[
        plot_df["status"] == status
    ]

    plt.scatter(
        subset["log2FoldChange"],
        subset["minus_log10_padj"],
        label=status,
        alpha=0.6,
        s=15
    )

# 6. ADD THRESHOLD LINES
plt.axvline(
    x=1,
    linestyle="--"
)

plt.axvline(
    x=-1,
    linestyle="--"
)

plt.axhline(
    y=-np.log10(0.05),
    linestyle="--"
)

# 7. LABEL EGFP
if "eGFP" in plot_df.index:

    egfp = plot_df.loc["eGFP"]

    plt.annotate(
        "eGFP",
        (
            egfp["log2FoldChange"],
            egfp["minus_log10_padj"]
        ),
        xytext=(5, 5),
        textcoords="offset points",
        fontsize=9
    )

# 8. LABEL AXES
plt.xlabel(
    "log2 Fold Change (mRNA-LNP vs saline)"
)

plt.ylabel(
    "-log10 adjusted p-value"
)

plt.title(
    "PyDESeq2 Differential Expression: "
    "mRNA-LNP vs Saline"
)

plt.legend()

plt.tight_layout()

# 9. SAVE
plt.savefig(
    "results/figures/pydeseq2_volcano.png",
    dpi=300
)

plt.close()

print("\nVolcano plot saved:")
print(
    "results/figures/pydeseq2_volcano.png"
)