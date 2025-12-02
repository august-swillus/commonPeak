#' Identify common ChIP peaks between biological conditions
#'
#' This function runs a 4-step pipeline:
#'   1. Find intersecting peaks across all replicates (bedtools multiinter)
#'   2. Count reads in those peaks 
#'   3. Normalize read counts 
#'   4. Identify commonly strong peaks between conditions (DESeq2)
#'
#' @param ... Named file paths to bed and bam files.
#' @param narrow Logical, TRUE for narrow peaks (default), FALSE for broad.
#' @param qvalue_cutoff Numeric, q-value cutoff for filtering peaks (default 0.05).
#' @param lfc_threshold Numeric, log2 fold change threshold for DESeq2 (default 1).
#' @param output_dir Character, directory for results (required).
#' @param paired Logical, TRUE if BAM files are paired-end; FALSE for single-end (default).
#' @param ma_plot create an MA plot (default)
#' @param results_csv create a results file listing the most common peaks (default)
#' @return A list containing DESeq2 results and the MA plot.
#' @export
commonPeak <- function(...,
                       narrow = TRUE,
                       qvalue_cutoff = 0.05,
                       lfc_threshold = 1,
                       output_dir = NULL,
                       paired = FALSE,
                       ma_plot = TRUE,
                       results_csv = TRUE) {

  args <- list(...)

  if (is.null(output_dir)) {
    stop("Please specify an output directory using output_dir = '/path/to/results/'.")
  }

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  bed_files <- args[grepl("bed", names(args))]
  bam_files <- args[grepl("bam", names(args))]

  if (length(bed_files) == 0 || length(bam_files) == 0) {
    stop("Please provide both bed and bam file paths.")
  }

# --- Parse and save manifests ---
    
  parse_group <- function(x) sub("^group_?([0-9]+)_.*", "group_\\1", x)
  bed_df <- data.frame(
    group = parse_group(names(bed_files)),
    sample = names(bed_files),
    path = unlist(bed_files),
    type = "peak",
    stringsAsFactors = FALSE
  )

  bam_df <- data.frame(
    group = parse_group(names(bam_files)),
    sample = tools::file_path_sans_ext(basename(unlist(bam_files))),
    path = unlist(bam_files),
    is_input = grepl("input", names(bam_files), ignore.case = TRUE),
    is_paired = paired,
    stringsAsFactors = FALSE
  )

  bed_manifest <- file.path(output_dir, "beds_manifest.tsv")
  bam_manifest <- file.path(output_dir, "bams_manifest.tsv")
  write.table(bed_df, bed_manifest, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(bam_df, bam_manifest, sep = "\t", quote = FALSE, row.names = FALSE)

  message("BED manifest written to: ", bed_manifest)
  message("BAM manifest written to: ", bam_manifest)

# --- Define script paths ---
    
  bash_script  <- "/data/cephfs-1/home/users/ausw10_c/work/CommonPeak/inst/scripts/step1_bedtools_multiinter.sh"
  first_R_script <- "/data/cephfs-1/home/users/ausw10_c/work/CommonPeak/inst/scripts/step2_count_reads.R"
  python_script  <- "/data/cephfs-1/home/users/ausw10_c/work/CommonPeak/inst/scripts/step3_normalize_counts.py"
  step4_script   <- "/data/cephfs-1/home/users/ausw10_c/work/CommonPeak/inst/scripts/step4_deseq2_analysis.R"

# --- Run pipeline ---
    
  message("Running Step 1: bedtools multiinter ...")
  system2("bash", c(bash_script, bed_manifest, output_dir))

  message("Running Step 2: read counting ...")
  peak_type <- ifelse(narrow, "narrow", "broad")
  paired_mode <- ifelse(paired, "paired", "single")
  system2("Rscript", c(first_R_script, bam_manifest, output_dir, peak_type, paired_mode))

  message("Running Step 3: normalization ...")
  system2("python3", c(python_script, bam_manifest, output_dir))

  message("Running Step 4: DESeq2 analysis ...")
  normalized_counts <- file.path(output_dir, "normalized_counts.csv")

  system2("Rscript", c(
    step4_script,
    normalized_counts,
    bam_manifest,
    output_dir,
    as.character(lfc_threshold),
    as.character(qvalue_cutoff),
    as.character(ma_plot),
    as.character(results_csv)
  ))

  message("Pipeline completed. Results written to: ", output_dir)
}



