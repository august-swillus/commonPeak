#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(Rsubread)
  library(GenomicAlignments)
  library(Rsamtools)
  library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript step2_count_reads.R <bam_manifest.tsv> <output_dir> <narrow|broad> [paired|single] [cores]")
}

bam_manifest <- args[1]
output_dir   <- args[2]
peak_type    <- args[3]
paired_mode  <- ifelse(length(args) >= 4, args[4], "single")
cores_arg    <- ifelse(length(args) >= 5, args[5], NA_character_)

if (!peak_type %in% c("narrow", "broad")) {
  stop("Error: Third argument must be either 'narrow' or 'broad'.")
}

# Determine paired-end status
is_paired <- paired_mode == "paired"
cat("Paired-end mode:", is_paired, "\n")

# Flank size setup
flank_size <- ifelse(peak_type == "narrow", 200, 500)
cat("Peak type:", peak_type, "- flank size set to ±", flank_size, "bp\n")

# Temporary directory
tmp_dir <- file.path(output_dir, "work", "tmp")
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
Sys.setenv(TMPDIR = tmp_dir)
options(tmpdir = tmp_dir)
cat("Temporary files will be written to:", tmp_dir, "\n")

# ---- Resolve cores (user -> SLURM -> detectCores -> fallback) ----
resolve_cores <- function(cores_arg, default = 4) {
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

num_cores <- resolve_cores(cores_arg, default = 4)
register(MulticoreParam(workers = num_cores))
cat("BiocParallel will use", num_cores, "core(s) (MulticoreParam).\n")

# Load Peaks
peak_file <- file.path(output_dir, "peak_in_all_samples.txt")
if (!file.exists(peak_file)) {
  stop(paste("Peak file not found in", output_dir, "- expected 'peak_in_all_samples.txt'"))
}

peak_data <- read.table(peak_file, header = FALSE, sep = "\t")
gr_peaks <- GRanges(seqnames = peak_data$V1, ranges = IRanges(start = peak_data$V2, end = peak_data$V3))
cat("Loaded", length(gr_peaks), "consensus peaks from:", peak_file, "\n")

# Load BAM manifest
bam_df <- read.table(bam_manifest, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
bam_files_all <- bam_df$path
names(bam_files_all) <- bam_df$sample
cat("Loaded", length(bam_files_all), "BAM files from:", bam_manifest, "\n")

# Identify summits using first BAM sample
cat("Identifying summits in first BAM sample:", names(bam_files_all)[1], "\n")
summit_bam_file <- bam_files_all[1]

find_summit_window_fast <- function(bam_file, peaks, flank_size) {
  cat("Processing", bam_file, "\n")

  aln <- readGAlignments(bam_file, use.names = TRUE)
  counts_per_base <- coverage(aln)
  peaks_with_cov <- peaks[seqnames(peaks) %in% names(counts_per_base)]

  summit_list <- BiocParallel::bplapply(seq_along(peaks_with_cov), function(i) {
    peak <- peaks_with_cov[i]
    chr <- as.character(seqnames(peak))
    if (!chr %in% names(counts_per_base)) return(NULL)

    cov_chr <- counts_per_base[[chr]]
    adj_start <- max(1, start(peak))
    adj_end   <- min(length(cov_chr), end(peak))
    peak_cov <- cov_chr[adj_start:adj_end]

    if (length(peak_cov) == 0 || all(peak_cov == 0)) return(NULL)

    max_pos <- which.max(as.vector(peak_cov))
    summit  <- adj_start + max_pos - 1
    summit_start <- max(summit - flank_size, 1)
    summit_end   <- summit + flank_size

    GRanges(seqnames = chr, ranges = IRanges(summit_start, summit_end))
  })

  Filter(Negate(is.null), summit_list)
}

summit_list <- find_summit_window_fast(summit_bam_file, gr_peaks, flank_size)
summit_gr <- do.call(c, summit_list)
cat("Summit-centered ±", flank_size, "bp windows created for", length(summit_gr), "peaks.\n")

# SAF conversion
saf <- data.frame(
  GeneID = paste0(seqnames(summit_gr), ":", start(summit_gr), "-", end(summit_gr)),
  Chr    = as.character(seqnames(summit_gr)),
  Start  = start(summit_gr),
  End    = end(summit_gr),
  Strand = "*",
  stringsAsFactors = FALSE
)

if (any(duplicated(saf$GeneID))) {
  dup_count <- sum(duplicated(saf$GeneID))
  cat("Warning:", dup_count, "duplicate peak IDs detected.\n")
  saf$GeneID <- make.unique(saf$GeneID)
}

# featureCounts
cat("Running featureCounts on", length(bam_files_all), "BAM files using", num_cores, "thread(s)...\n")

fc_results <- featureCounts(
  files = bam_files_all,
  annot.ext = saf,
  isGTFAnnotationFile = FALSE,
  useMetaFeatures = FALSE,
  allowMultiOverlap = TRUE,
  isPairedEnd = is_paired,
  strandSpecific = 0,
  nthreads = num_cores,
  tmpDir = tmp_dir,
  verbose = FALSE
)

# Save results
count_matrix <- fc_results$counts
colnames(count_matrix) <- names(bam_files_all)
rownames(count_matrix) <- saf$GeneID

count_file <- file.path(output_dir, "read_counts_matrix.csv")
write.csv(count_matrix, count_file)
cat("Count matrix saved to:", count_file, "\n")
