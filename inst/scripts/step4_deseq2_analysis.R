#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
  stop("Usage: Rscript step4_deseq2_analysis.R <normalized_counts.csv> <bam_manifest.tsv> <output_dir> <lfcThreshold> <alpha> <ma_plot> <results_csv> [cores]")
}

count_file    <- args[1]
bam_manifest  <- args[2]
output_dir    <- args[3]
lfc_threshold <- as.numeric(args[4])
alpha_cutoff  <- as.numeric(args[5])
ma_plot       <- as.logical(args[6])
results_csv   <- as.logical(args[7])
cores_arg     <- ifelse(length(args) >= 8, args[8], NA_character_)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("Running DESeq2 with lfcThreshold =", lfc_threshold,
    "and q-value cutoff =", alpha_cutoff, "\n")
cat("MA plot:", ma_plot, "| Results CSV:", results_csv, "\n")

# ---- Resolve cores (user -> SLURM -> detectCores -> fallback) ----
resolve_cores <- function(cores_arg, default = 8) {
  if (!is.na(cores_arg) && nzchar(cores_arg)) {
    cval <- suppressWarnings(as.integer(cores_arg))
    if (!is.na(cval) && cval >= 1) return(cval)
    stop("Invalid cores argument. Provide a single integer >= 1.")
  }

  slurm <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")))
  if (!is.na(slurm) && slurm >= 1) return(slurm)

  dc <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (!is.na(dc) && dc >= 1) return(dc)

  default
}

num_cores <- resolve_cores(cores_arg, default = 8)
register(MulticoreParam(workers = num_cores))
cat("Using", num_cores, "core(s) for DESeq2 (MulticoreParam).\n")

# Load input data
cts <- read.csv(count_file, row.names = 1, check.names = FALSE)
cts <- round(cts)
bam_df <- read.table(bam_manifest, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Verify sample overlap
common_samples <- intersect(colnames(cts), bam_df$sample)
if (length(common_samples) < 2) {
  stop("Error: Fewer than 2 overlapping sample names between normalized counts and BAM manifest.")
}
cts <- cts[, common_samples, drop = FALSE]
bam_df <- bam_df[bam_df$sample %in% common_samples, ]

# Define groups
groups <- unique(bam_df$group)
if (length(groups) != 2) {
  stop(paste("Exactly two groups required for DESeq2 comparison. Found:", paste(groups, collapse = ", ")))
}

coldata <- data.frame(condition = factor(bam_df$group, levels = groups))
rownames(coldata) <- bam_df$sample
cat("Comparing groups:", paste(groups, collapse = " vs "), "\n")

# Prepare DESeq2 dataset
keep <- rowSums(cts) >= 10
cts_filtered <- cts[keep, , drop = FALSE]

dds <- DESeqDataSetFromMatrix(
  countData = cts_filtered,
  colData   = coldata,
  design    = ~ condition
)

# Run DESeq2
dds <- DESeq(dds, parallel = TRUE)
contrast <- c("condition", groups[2], groups[1])

res <- results(
  dds,
  contrast      = contrast,
  lfcThreshold  = lfc_threshold,
  altHypothesis = "lessAbs",
  parallel      = TRUE,
  alpha         = alpha_cutoff
)

resLFC <- lfcShrink(dds, coef = 2, type = "apeglm", res = res)

cat("\nDESeq2 analysis complete.\n")
cat("Comparison:", groups[2], "vs", groups[1], "\n\n")

# Print summary
total_nonzero <- sum(resLFC$baseMean > 0, na.rm = TRUE)
up <- sum(resLFC$log2FoldChange > 0 &
            resLFC$log2FoldChange < lfc_threshold &
            resLFC$padj < alpha_cutoff, na.rm = TRUE)
down <- sum(resLFC$log2FoldChange < 0 &
              resLFC$log2FoldChange > -lfc_threshold &
              resLFC$padj < alpha_cutoff, na.rm = TRUE)

cat("out of", total_nonzero, "with nonzero total read count\n")
cat("Fold change > 0 < ", lfc_threshold, ": ", up, "\n", sep = "")
cat("Fold change < 0 > -", lfc_threshold, ": ", down, "\n", sep = "")

# Order and filter results
resOrdered <- resLFC[order(resLFC$pvalue), ]
resSig <- subset(resOrdered, padj < alpha_cutoff)

# Add genomic coordinates to results
if (nrow(resSig) > 0) {
  resSig$peak_id <- rownames(resSig)
  coords <- do.call(rbind, strsplit(resSig$peak_id, "[:-]"))
  colnames(coords) <- c("chr", "start", "end")
  resSig <- cbind(coords, resSig)
  rownames(resSig) <- NULL
}

if (results_csv) {
  out_csv <- file.path(output_dir, "DESeq2_results_with_coords.csv")
  write.csv(resSig, out_csv, row.names = FALSE)
  cat("\nFiltered DESeq2 results saved to:", out_csv, "\n")
} else {
  cat("\nSkipping CSV export (results_csv = FALSE)\n")
}

# Generate MA plot
if (ma_plot) {
  ylim <- c(-3, 3)
  ma_plot_file <- file.path(output_dir, "DESeq2_MA_plot.png")
  png(ma_plot_file, width = 800, height = 600, res = 120)
  plotMA(resLFC, ylim = ylim)
  abline(h = c(-lfc_threshold, lfc_threshold), col = "dodgerblue", lwd = 2)
  dev.off()
  cat("MA plot saved to:", ma_plot_file, "\n")
} else {
  cat("Skipping MA plot generation (ma_plot = FALSE)\n")
}

cat("\nStep 4 complete.\n")
