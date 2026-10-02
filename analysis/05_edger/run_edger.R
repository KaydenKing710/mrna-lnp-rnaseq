library(edgeR)

# 1 LOAD THE SAME FILTERED COUNTS
# USED BY PyDESeq2
counts <- read.csv(
  "data/processed/filtered_counts.csv",
  row.names = 1,
  check.names = FALSE
)

cat("Input matrix:\n")
cat(
  nrow(counts),
  "samples x",
  ncol(counts),
  "genes\n"
)

# 2 TRANSPOSE FOR edgeR
# PyDESeq2 uses:
# samples x genes
#
# edgeR uses:
# genes x samples

counts <- t(counts)

cat("\nMatrix for edgeR:\n")
cat(
  nrow(counts),
  "genes x",
  ncol(counts),
  "samples\n"
)

# 3 CHECK SAMPLE ORDER
expected_samples <- c(
  "saline1",
  "saline2",
  "saline3",
  "saline4",
  "saline5",
  "mRNA-eGFP1",
  "mRNA-eGFP2",
  "mRNA-eGFP3",
  "mRNA-eGFP4",
  "mRNA-eGFP5"
)

cat("\nSample names:\n")
print(colnames(counts))

stopifnot(
  identical(
    colnames(counts),
    expected_samples
  )
)

cat("\nSample order check: PASSED\n")

# 4 DEFINE EXPERIMENTAL GROUPS
group <- factor(
  c(
    rep("control", 5),
    rep("mRNA_LNP", 5)
  ),
  levels = c(
    "control",
    "mRNA_LNP"
  )
)

cat("\nExperimental groups:\n")
print(group)

# 5 CREATE DGEList
y <- DGEList(
  counts = counts,
  group = group
)


# 6 TMM NORMALIZATION
y <- calcNormFactors(
  y,
  method = "TMM"
)

cat("\nLibrary information after TMM normalization:\n")
print(y$samples)

# 7 CREATE DESIGN MATRIX
design <- model.matrix(
  ~ group
)

cat("\nDesign matrix:\n")
print(design)

# 8 ESTIMATE DISPERSION
y <- estimateDisp(
  y,
  design
)

cat("\nDispersion estimation completed.\n")

# 9 FIT QUASI-LIKELIHOOD MODEL
fit <- glmQLFit(
  y,
  design
)

cat(
  "\nQuasi-likelihood model fitting completed.\n"
)

# 10 TEST mRNA_LNP VS CONTROL
qlf <- glmQLFTest(
  fit,
  coef = "groupmRNA_LNP"
)

# 11. EXTRACT ALL RESULTS
results <- topTags(
  qlf,
  n = Inf,
  sort.by = "none"
)$table

cat("\nFirst 5 edgeR results:\n")
print(
  head(results, 5)
)

# 12. SAVE ALL RESULTS
write.csv(
  results,
  "results/tables/edger_all_results.csv",
  row.names = TRUE
)

# 13 DEFINE SIGNIFICANT DEGs
significant <- results[
  results$FDR < 0.05 &
  abs(results$logFC) >= 1,
]

upregulated <- significant[
  significant$logFC >= 1,
]

downregulated <- significant[
  significant$logFC <= -1,
]

# 14. PRINT DEG SUMMARY
cat("\n")
cat("====================================\n")
cat("edgeR DEG SUMMARY\n")
cat("====================================\n")

cat(
  "Genes tested:",
  nrow(results),
  "\n"
)

cat(
  "Significant DEGs:",
  nrow(significant),
  "\n"
)

cat(
  "Upregulated:",
  nrow(upregulated),
  "\n"
)

cat(
  "Downregulated:",
  nrow(downregulated),
  "\n"
)
# 15. CHECK eGFP
if ("eGFP" %in% rownames(results)) {

  cat("\neGFP edgeR result:\n")

  print(
    results["eGFP", ]
  )
}

# 16. SAVE SIGNIFICANT DEGs
write.csv(
  significant,
  "results/tables/edger_significant_DEGs.csv",
  row.names = TRUE
)

cat(
  "\nedgeR results saved successfully.\n"
)