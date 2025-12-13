#' Identify common ChIP peaks between biological conditions
#'
#' This function runs a 4-step pipeline:
#'   1. Find intersecting peaks across all replicates (bedtools multiinter)
#'   2. Count reads in those peaks
#'   3. Normalize read counts
#'   4. Identify commonly strong peaks between conditions (DESeq2 equivalence test)
#'
#' @param ... Named file paths to bed and bam files.
#' @param narrow Logical, TRUE for narrow peaks (default), FALSE for broad.
#' @param qvalue_cutoff Numeric, q-value cutoff for filtering peaks (default 0.05).
#' @param lfc_threshold Numeric, log2 fold change threshold for DESeq2 (default 1).
#' @param output_dir Character, directory for results (required).
#' @param paired Logical, TRUE if BAM files are paired-end; FALSE for single-end (default).
#' @param ma_plot Logical, create an MA plot (default TRUE).
#' @param results_csv Logical, create a results CSV file (default TRUE).
#' @param cores Integer or NULL. Number of CPU cores to use for parallel steps (Steps 2 and 4).
#'   If NULL, scripts will try SLURM_CPUS_PER_TASK, then parallel::detectCores(), then fall back.
#' @return Invisibly returns NULL (pipeline writes results to output_dir).
#' @export
commonPeak <- function(...,
                       narrow = TRUE,
                       qvalue_cutoff = 0.05,
                       lfc_threshold = 1,
                       output_dir = NULL,
                       paired = FALSE,
                       ma_plot = TRUE,
                       results_csv = TRUE,
                       cores = NULL) {

  args <- list(...)

  if (is.null(output_dir) || !nzchar(output_dir)) {
    stop("Please specify an output directory using output_dir = '/path/to/results/'.")
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  if (!is.null(cores)) {
    if (!is.numeric(cores) || length(cores) != 1 || is.na(cores) || cores < 1) {
      stop("cores must be a single integer >= 1 (or NULL).")
    }
    cores <- as.integer(cores)
  }

  # -------------------------
  # Logging setup
  # -------------------------
  log_file <- file.path(output_dir, "commonPeak_pipeline.log")

  timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_write <- function(lines, append = TRUE) {
    if (length(lines) == 0) return(invisible(NULL))
    writeLines(lines, con = log_file, sep = "\n", useBytes = TRUE)
    if (append) writeLines("", con = log_file, sep = "\n", useBytes = TRUE)
    invisible(NULL)
  }
  log_line <- function(...) {
    msg <- paste0("[", timestamp(), "] ", paste0(..., collapse = ""))
    writeLines(msg, con = log_file, sep = "\n", useBytes = TRUE)
    invisible(NULL)
  }
  log_block <- function(header, lines) {
    log_line(header)
    if (length(lines)) {
      writeLines(lines, con = log_file, sep = "\n", useBytes = TRUE)
    }
    writeLines("", con = log_file, sep = "\n", useBytes = TRUE)
    invisible(NULL)
  }

  # Create/overwrite log file with header
  writeLines(character(0), con = log_file)
  log_line("commonPeak pipeline log started")
  log_line("output_dir: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE))

  # Record basic environment/package info (best-effort)
  pkg_version <- tryCatch(as.character(utils::packageVersion("commonPeak")), error = function(e) NA_character_)
  log_line("commonPeak version: ", pkg_version)
  log_line("R version: ", R.version.string)
  log_line("Platform: ", R.version$platform)
  log_line("Working directory: ", getwd())

  # Record user parameters
  log_block("Parameters:", c(
    paste0("  narrow = ", narrow),
    paste0("  qvalue_cutoff = ", qvalue_cutoff),
    paste0("  lfc_threshold = ", lfc_threshold),
    paste0("  paired = ", paired),
    paste0("  ma_plot = ", ma_plot),
    paste0("  results_csv = ", results_csv),
    paste0("  cores = ", ifelse(is.null(cores), "NULL", as.character(cores)))
  ))

  # Record provided arguments (names + paths)
  if (length(args)) {
    arg_lines <- unlist(lapply(names(args), function(nm) paste0("  ", nm, " = ", args[[nm]])))
    log_block("Input arguments (...):", arg_lines)
  }

  # Identify bed/bam inputs by argument names
  bed_files <- args[grepl("bed", names(args), ignore.case = TRUE)]
  bam_files <- args[grepl("bam", names(args), ignore.case = TRUE)]

  if (length(bed_files) == 0 || length(bam_files) == 0) {
    log_line("ERROR: Missing BED and/or BAM inputs (argument names must contain 'bed' / 'bam').")
    stop("Please provide both bed and bam file paths (argument names must contain 'bed' / 'bam'). ",
         "See log: ", log_file)
  }

  # --- Parse and save manifests ---
  parse_group <- function(x) sub("^group_?([0-9]+)_.*", "group_\\1", x)

  bed_df <- data.frame(
    group  = parse_group(names(bed_files)),
    sample = names(bed_files),
    path   = unlist(bed_files),
    type   = "peak",
    stringsAsFactors = FALSE
  )

  bam_df <- data.frame(
    group     = parse_group(names(bam_files)),
    sample    = tools::file_path_sans_ext(basename(unlist(bam_files))),
    path      = unlist(bam_files),
    is_input  = grepl("input", names(bam_files), ignore.case = TRUE),
    is_paired = paired,
    stringsAsFactors = FALSE
  )

  bed_manifest <- file.path(output_dir, "beds_manifest.tsv")
  bam_manifest <- file.path(output_dir, "bams_manifest.tsv")
  utils::write.table(bed_df, bed_manifest, sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(bam_df, bam_manifest, sep = "\t", quote = FALSE, row.names = FALSE)

  message("BED manifest written to: ", bed_manifest)
  message("BAM manifest written to: ", bam_manifest)
  log_line("BED manifest written to: ", bed_manifest)
  log_line("BAM manifest written to: ", bam_manifest)

  # Log manifest contents
  log_block("beds_manifest.tsv contents:", readLines(bed_manifest, warn = FALSE))
  log_block("bams_manifest.tsv contents:", readLines(bam_manifest, warn = FALSE))

  # --- Locate scripts inside the installed package ---
  bash_script   <- system.file("scripts", "step1_bedtools_multiinter.sh", package = "commonPeak")
  step2_script  <- system.file("scripts", "step2_count_reads.R",          package = "commonPeak")
  python_script <- system.file("scripts", "step3_normalize_counts.py",    package = "commonPeak")
  step4_script  <- system.file("scripts", "step4_deseq2_analysis.R",      package = "commonPeak")

  log_block("Resolved script paths:", c(
    paste0("  step1: ", bash_script),
    paste0("  step2: ", step2_script),
    paste0("  step3: ", python_script),
    paste0("  step4: ", step4_script)
  ))

  if (bash_script == "" || step2_script == "" || python_script == "" || step4_script == "") {
    log_line("ERROR: One or more pipeline scripts not found via system.file().")
    stop("One or more pipeline scripts were not found via system.file(). ",
         "Ensure they exist under inst/scripts/ and are included in the package build. ",
         "See log: ", log_file)
  }

  # Optional: record tool versions (best-effort; do not fail pipeline if unavailable)
  run_version_probe <- function(cmd, cmd_args, label) {
    out <- tryCatch(system2(cmd, cmd_args, stdout = TRUE, stderr = TRUE), error = function(e) NULL)
    if (is.null(out)) {
      log_line(label, ": version probe failed (command could not be executed).")
    } else {
      log_block(paste0(label, " version probe output:"), out)
    }
  }
  run_version_probe("bedtools", c("--version"), "bedtools")
  run_version_probe("python3", c("--version"), "python3")
  run_version_probe("python3", c("-c", "import pandas, numpy; print('pandas', pandas.__version__); print('numpy', numpy.__version__)"),
                    "python3 (pandas/numpy)")

  # --- Helper: run commands and stop on failure (also logs stdout/stderr) ---
  run_cmd <- function(command, cmd_args, step_name) {
    log_line("----- START ", step_name, " -----")
    log_line("Command: ", command, " ", paste(cmd_args, collapse = " "))

    out <- system2(command, cmd_args, stdout = TRUE, stderr = TRUE)
    status <- attr(out, "status")
    exit_code <- if (is.null(status)) 0L else as.integer(status)

    # Write full output to log
    log_block(paste0(step_name, " output:"), out)

    # Also print to console so user sees progress as before
    if (length(out)) cat(paste0(out, collapse = "\n"), "\n")

    if (!identical(exit_code, 0L)) {
      log_line("----- FAIL ", step_name, " (exit code ", exit_code, ") -----")
      stop(step_name, " failed (exit code ", exit_code, "). See log file: ", log_file)
    }

    log_line("----- END ", step_name, " (exit code 0) -----")
    invisible(out)
  }

  # -------------------------
  # Run pipeline
  # -------------------------
  message("Log file: ", log_file)
  log_line("Log file: ", log_file)

  message("Running Step 1: bedtools multiinter ...")
  run_cmd("bash", c(bash_script, bed_manifest, output_dir), "Step 1 (bedtools multiinter)")

  message("Running Step 2: read counting ...")
  peak_type   <- ifelse(isTRUE(narrow), "narrow", "broad")
  paired_mode <- ifelse(isTRUE(paired), "paired", "single")

  step2_args <- c(step2_script, bam_manifest, output_dir, peak_type, paired_mode)
  if (!is.null(cores)) step2_args <- c(step2_args, as.character(cores))
  run_cmd("Rscript", step2_args, "Step 2 (read counting)")

  message("Running Step 3: normalization ...")
  run_cmd("python3", c(python_script, bam_manifest, output_dir), "Step 3 (normalization)")

  message("Running Step 4: DESeq2 analysis ...")
  normalized_counts <- file.path(output_dir, "normalized_counts.csv")

  step4_args <- c(
    step4_script,
    normalized_counts,
    bam_manifest,
    output_dir,
    as.character(lfc_threshold),
    as.character(qvalue_cutoff),
    as.character(ma_plot),
    as.character(results_csv)
  )
  if (!is.null(cores)) step4_args <- c(step4_args, as.character(cores))
  run_cmd("Rscript", step4_args, "Step 4 (DESeq2 analysis)")

  message("Pipeline completed. Results written to: ", output_dir)
  log_line("Pipeline completed successfully. Results written to: ", output_dir)
  log_line("commonPeak pipeline log finished")

  invisible(NULL)
}
