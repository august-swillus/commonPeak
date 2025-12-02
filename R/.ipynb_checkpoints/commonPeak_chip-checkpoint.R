#' Identify common ChIP peaks between biological conditions
#'
#' This function runs a 4-step pipeline:
#'   1. Find intersecting peaks across all replicates (bedtools multiinter)
#'   2. Count reads in those peaks 
#'   3. Normalize read counts 
#'   4. Identify commonly strong peaks between conditions (DESeq2)
#'
#' @param ... Named file paths to bed and bam files.
#'   Example:
#'   \code{commonPeak(group1_bed_s1="path/to/file1.bed", group1_bed_s2="path/to/file2.bed", group2_bed_s1="path/to/file3.bed", group1_bam_s1="path/to/file1.bam", group1_bam_s2="path/to/file2.bam", group2_bam_s1="path/to/file3.bam", narrow=TRUE)}
#' @param narrow Logical, TRUE for narrow peaks (default), FALSE for broad.
#' @param qvalue_cutoff Numeric, q-value cutoff for filtering peaks (default 0.05).
#' @param output_dir Character, path to directory where results should be saved (required).
#' @return A list containing DESeq2 results and the MA plot.
#' @export
#' @examples
#' \dontrun{
#' res <- commonPeak(
#'  group1_bed_s1 = "path/to/ChIP_rep1.bed",
#'  group1_bed_s2 = "path/to/ChIP_rep2.bed",
#'  group1_bam_s1 = "path/to/ChIP_rep1.bam",
#'  group1_bam_s2 = "path/to/ChIP_rep2.bam",
#'  group1_input_bam = "path/to/input_group1.bam",   
#'  group2_bed_s1 = "path/to/iChIP_rep1.bed",
#'  group2_bed_s2 = "path/to/iChIP_rep2.bed",
#'  group2_bam_s1 = "path/to/iChIP_rep1.bam",
#'  group2_bam_s2 = "path/to/iChIP_rep2.bam",
#'  group2_input_bam = "path/to/input_group2.bam",   
#'  output_dir = "path/to/results"
#')

#' }

commonPeak <- function(...,
                       narrow = TRUE,
                       qvalue_cutoff = 0.05,
                       lfc_threshold = 0.75,
                       output_dir = NULL) {
  
  args <- list(...)
  
  # Validate output_dir
  if (is.null(output_dir)) {
    stop("Please specify an output directory using the 'output_dir' argument, e.g., output_dir = '/path/to/results/'.")
  }
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Separate bed and bam files
  bed_files <- args[grepl("bed", names(args))]
  bam_files <- args[grepl("bam", names(args))]
  
  if (length(bed_files) == 0 || length(bam_files) == 0) {
    stop("Please provide both bed and bam file paths.")
  }
  
  # Parse group info
  parse_group <- function(x) sub("^group_?([0-9]+)_.*", "group_\\1", x)
  bed_groups <- parse_group(names(bed_files))
  bam_groups <- parse_group(names(bam_files))
  
  bed_df <- data.frame(
    group = bed_groups,
    sample = names(bed_files),
    path = unlist(bed_files),
    type = "peak",
    stringsAsFactors = FALSE
  )
  
  # --- Identify BAMs ---
  bam_names <- names(bam_files)

  is_input <- grepl("input", bam_names, ignore.case = TRUE)
  groups <- parse_group(bam_names)

  bam_df <- data.frame(
  group = bam_groups,
  sample = tools::file_path_sans_ext(basename(unlist(bam_files))),  
  path = unlist(bam_files),
  is_input = grepl("input", names(bam_files), ignore.case = TRUE),
  stringsAsFactors = FALSE
)

  # --- Validation ---
  if (any(is_input)) {
    message("Detected input BAMs: ", paste(bam_df$sample[bam_df$is_input], collapse = ", "))
  } else {
    warning(" No input BAMs detected. Normalization will proceed without subtraction.")
  }

  
  # Write manifests
  bed_manifest <- file.path(output_dir, "beds_manifest.tsv")
  bam_manifest <- file.path(output_dir, "bams_manifest.tsv")
  write.table(bed_df, bed_manifest, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(bam_df, bam_manifest, sep = "\t", quote = FALSE, row.names = FALSE)
  
  message("BED manifest written to: ", bed_manifest)
  message("BAM manifest written to: ", bam_manifest)
  
  # Define script paths
  bash_script  <- "/data/cephfs-1/home/users/ausw10_c/work/CommonPeak/inst/scripts/step1_bedtools_multiinter.sh"
  first_R_script <- "/data/cephfs-1/home/users/ausw10_c/work/CommonPeak/inst/scripts/step2_count_reads.R"
  python_script  <- "/data/cephfs-1/home/users/ausw10_c/work/CommonPeak/inst/scripts/step3_normalize_counts.py"
  step4_script   <- "/data/cephfs-1/home/users/ausw10_c/work/CommonPeak/inst/scripts/step4_deseq2_analysis.R"
  
  # Step 1
  message("Running Step 1: bedtools multiinter ...")
  system2("bash", c(bash_script, bed_manifest, output_dir))
  
  # Step 2
  message("Running Step 2: read counting ...")
  peak_type <- ifelse(narrow, "narrow", "broad")
  system2("Rscript", c(first_R_script, bam_manifest, output_dir, peak_type))
  
  # Step 3
  message("Running Step 3: normalization ...")
  system2("python3", c(python_script, bam_manifest, output_dir))
  
  # Step 4
  message("Running Step 4: DESeq2 analysis ...")
  normalized_counts <- file.path(output_dir, "normalized_counts.csv") # adjust if needed
  
  system2("Rscript", c(
    step4_script,
    normalized_counts,
    bam_manifest,
    output_dir,
    as.character(lfc_threshold),
    as.character(qvalue_cutoff)
  ))
  
  message("Pipeline completed. Results written to: ", output_dir)
}


#test if commonpeak function correctly calls multiinter script leading to output files being written 
# ---- Test run of commonPeak function ----
res <- commonPeak(
  group1_bed_s1 = "/data/cephfs-1/home/users/ausw10_c/work/ichip_vs_chip/results/MACS3_blacklists_removed/k27ac-30M/k27ac-30MChIP_MCAS_H3K27Aca_S13.blacklists_removed_30M_peaks.narrowPeak",
  group1_bed_s2 = "/data/cephfs-1/home/users/ausw10_c/work/ichip_vs_chip/results/MACS3_blacklists_removed/k27ac-30M/k27ac-30MChIP_MCAS_H3K27Acb_S14.blacklists_removed_30M_peaks.narrowPeak", 
  group2_bed_s1 = "/data/cephfs-1/home/users/ausw10_c/work/ichip_vs_chip/results/MACS3_blacklists_removed/k27ac-30M/k27ac-30MiChIP_Fresh_MCAS_H3K27Aca_S73.blacklists_removed_30M_peaks.narrowPeak",
  group2_bed_s2 = "/data/cephfs-1/home/users/ausw10_c/work/ichip_vs_chip/results/MACS3_blacklists_removed/k27ac-30M/k27ac-30MiChIP_Fresh_MCAS_H3K27Acb_S74.blacklists_removed_30M_peaks.narrowPeak",
    
  group1_bam_s1 = "/data/cephfs-1/home/users/ausw10_c/work/ichip_vs_chip/results/picard_blacklists_removed/picard_deduplicated_bams/downsampled_mixed/ChIP_MCAS_H3K27Aca_S13.blacklists_duplicates_actually_removed_barcode_aware_30M.bam",
  group1_bam_s2 = "/data/cephfs-1/home/users/ausw10_c/work/ichip_vs_chip/results/picard_blacklists_removed/picard_deduplicated_bams/downsampled_mixed/ChIP_MCAS_H3K27Acb_S14.blacklists_duplicates_actually_removed_barcode_aware_30M.bam",
  group2_bam_s1 = "/data/cephfs-1/home/users/ausw10_c/work/ichip_vs_chip/results/picard_blacklists_removed/picard_deduplicated_bams/downsampled_mixed/iChIP_Fresh_MCAS_H3K27Aca_S73.blacklists_duplicates_actually_removed_barcode_aware_30M.bam",
  group2_bam_s2 = "/data/cephfs-1/home/users/ausw10_c/work/ichip_vs_chip/results/picard_blacklists_removed/picard_deduplicated_bams/downsampled_mixed/iChIP_Fresh_MCAS_H3K27Acb_S74.blacklists_duplicates_actually_removed_barcode_aware_30M.bam",

  group1_input_bam = "/data/cephfs-1/home/users/ausw10_c/work/ichip_vs_chip/results/picard_blacklists_removed/picard_deduplicated_bams/downsampled_mixed/ChIP_MCAS_ChIP_input_50ng_S25.blacklists_duplicates_actually_removed_barcode_aware_30M.bam",
  group2_input_bam = "/data/cephfs-1/home/users/ausw10_c/work/ichip_vs_chip/results/picard_blacklists_removed/picard_deduplicated_bams/downsampled_mixed/iChIP_Fresh_MCAS_Input_DNA_25ng_S27.blacklists_duplicates_actually_removed_barcode_aware_30M.bam",
    
  output_dir = "/data/cephfs-1/home/users/ausw10_c/work/CommonPeak/test_output/"
)


