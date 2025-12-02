#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(BiocParallel)
})

# ============================
# Parse Command-Line Arguments
# ============================

count_file    <- "/data/cephfs-1/home/users/ausw10_c/work/CommonPeak/test_output_mcf7/normalized_counts.csv"
bam_manifest  <- "/data/cephfs-1/home/users/ausw10_c/work/CommonPeak/test_output_mcf7/bams_manifest.tsv"
output_dir    <- "/data/cephfs-1/home/users/ausw10_c/work/CommonPeak/test_output_mcf7"
lfc_threshold <- 1
alpha_cutoff  <- 0.05


dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("Running DESeq2 with lfcThreshold =", lfc_threshold,
    "and p-value cutoff (alpha) =", alpha_cutoff, "\n")

# ============================
# Parallel setup
# ============================
num_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 8))
register(MulticoreParam(num_cores))
cat("Using", num_cores, "cores for DESeq2.\n")

# ============================
# Load input data
# ============================
cts <- read.csv(count_file, row.names = 1, check.names = FALSE)
cts <- round(cts)
bam_df <- read.table(bam_manifest, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Verify sample overlap and align orders EXACTLY to counts
samples <- colnames(cts)
if (!all(samples %in% bam_df$sample)) {
  stop("Error: Some count matrix columns are not present in BAM manifest 'sample' column.")
}
bam_df <- bam_df[match(samples, bam_df$sample), , drop = FALSE]
stopifnot(identical(bam_df$sample, samples))

# ============================
# Define groups (factor levels set BEFORE creating dds)
# ============================
groups <- unique(bam_df$group)
if (length(groups) != 2) {
  stop(paste("Exactly two groups required for DESeq2 comparison. Found:", paste(groups, collapse = ", ")))
}
coldata <- data.frame(condition = factor(bam_df$group, levels = groups),
                      row.names = bam_df$sample)
cat("Comparing groups (factor levels):", paste(levels(coldata$condition), collapse = " vs "), "\n")

# ============================
# Filter and create DESeq dataset
# ============================
keep <- rowSums(cts) >= 10
cts_filtered <- cts[keep, , drop = FALSE]

dds <- DESeqDataSetFromMatrix(
  countData = cts_filtered,
  colData   = coldata,
  design    = ~ condition
)

# ============================
# Run DESeq2
# ============================
dds <- DESeq(dds, parallel = TRUE)

# Results: same options as your working script
res <- results(
  dds,
  lfcThreshold   = lfc_threshold,
  altHypothesis  = "lessAbs",
  parallel       = TRUE,
  alpha          = alpha_cutoff
)

# Pick the correct coef for lfcShrink automatically
rn <- resultsNames(dds)
# Expect something like: "Intercept", "condition_group2_vs_group1"
coef_idx <- grep("^condition_.*_vs_.*$", rn)
if (length(coef_idx) != 1) {
  stop("Could not uniquely determine contrast coefficient for lfcShrink(). resultsNames(dds):\n",
       paste(rn, collapse = "\n"))
}

resLFC <- lfcShrink(dds, coef = coef_idx, type = "apeglm", res = res)

cat("\nDESeq2 analysis complete.\n")

# ============================
# Summary (exactly like your simple script)
# ============================
summary(resLFC)

# ============================
# Order by p-value and filter by padj < alpha
# ============================
resOrdered <- resLFC[order(resLFC$pvalue), ]
resSig <- subset(resOrdered, padj < alpha_cutoff)

# ============================
# Add genomic coordinates (robust to 0 or 1 rows)
# ============================
if (nrow(resSig) > 0) {
  resSig$peak_id <- rownames(resSig)
  coords_list <- strsplit(resSig$peak_id, "[:-]")
  coords <- do.call(rbind, lapply(coords_list, function(x) {
    if (length(x) == 3) x else c(NA, NA, NA)
  }))
  coords <- as.data.frame(coords, stringsAsFactors = FALSE)
  colnames(coords) <- c("chr", "start", "end")
  resSig <- cbind(coords, resSig)
  rownames(resSig) <- NULL
} else {
  warning("No significant peaks found — skipping coordinate parsing.")
}

# ============================
# Save results (as in your example)
# ============================
out_csv <- file.path(output_dir, "DESeq2_results_with_coords.csv")
write.csv(as.data.frame(resSig), out_csv, row.names = FALSE)
cat("\n Results written to:", out_csv, "\n")

# Save R objects for reproducibility
saveRDS(resLFC,file.path(output_dir, "deseq2_results.rds"))

# ============================
# MA Plot (exactly like your example; threshold lines use user value)
# ============================
ylim <- c(-3, 3)
ma_plot_file <- file.path(output_dir, "DESeq2_MA_plot.png")
png(ma_plot_file, width = 800, height = 600, res = 120)
plotMA(resLFC, ylim = ylim)
abline(h = c(-lfc_threshold, lfc_threshold), col = "dodgerblue", lwd = 2)
dev.off()
cat(" MA plot saved to:", ma_plot_file, "\n")

cat("\nStep 4 complete.\n")
