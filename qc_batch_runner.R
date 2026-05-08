#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# qc_batch_runner.R
#
# Render sample_QC.Rmd for every sample listed in a sample sheet.
#
# USAGE (interactive or on a compute node):
#   Rscript qc_batch_runner.R --sample_sheet samples.csv --output_dir results/QC
#
# SAMPLE SHEET FORMAT (CSV, with header):
#   sample_name,cellranger_dir,batch
#   Control_1,/path/to/cellranger/count-Control-1/outs,A
#   Treatment_1,/path/to/cellranger/count-Treatment-1/outs,A
#
# The script expects CellRanger's standard output layout, i.e.:
#   <cellranger_dir>/outs/filtered_feature_bc_matrix/
#   <cellranger_dir>/outs/raw_feature_bc_matrix/
# If your paths already point to the `outs` subfolder, set --outs_subdir ""
#
# Directory names are auto-detected when they differ slightly from standard
# names, as long as they contain a valid 10x MEX matrix layout.
#
# OPTIONAL ARGUMENTS:
#   --template      Path to the .Rmd template
#                   (default: same directory as this script)
#   --merge_template Path to merge_analysis.Rmd template
#                   (default: same directory as this script)
#   --outs_subdir   Subdirectory appended to cellranger_dir before
#                   filtered/raw paths (default: "outs")
#   --ncores        Number of samples to render in parallel (default: 1)
#   --merge_only    If TRUE, skip rendering and only merge existing
#                   <sample_name>_QC.rds files (default: FALSE)
#   --skip_merge    If TRUE, skip automatic merge after successful
#                   multi-sample rendering (default: FALSE)
#   --run_integration If TRUE, render integrate_RNA.Rmd after a successful
#                   merge (default: FALSE)
#   --integration_template Path to integrate_RNA.Rmd template
#                   (default: same directory as this script)
#   --integration_only If TRUE, run only integration on an existing
#                   merged_QC.rds (default: FALSE)
#   --integration_level Required when integration runs: Batch or Sample
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(rmarkdown)
  library(parallel)
})

# ── CLI argument parsing (base R; avoids external parser dependency) ──────────
print_usage_and_exit <- function(status = 0L) {
  cat(
    "Usage:\n",
    "  Rscript qc_batch_runner.R --sample_sheet <file.csv> [options]\n\n",
    "Options:\n",
    "  --sample_sheet <path>   Path to CSV sample sheet [required]\n",
    "  --output_dir <path>     Directory for HTML reports and RDS files [default: QC]\n",
    "  --template <path>       Path to sample_QC.Rmd [default: same dir as script]\n",
    "  --merge_template <path> Path to merge_analysis.Rmd [default: same dir as script]\n",
    "  --outs_subdir <name>    Subdirectory inside cellranger_dir [default: outs]\n",
    "  --ncores <int>          Number of samples to render in parallel [default: 1]\n",
    "  --sample_name <name>    Process only one sample_name from sheet\n",
    "  --merge_only <TRUE/FALSE>  If TRUE, skip individual sample QC and only merge existing RDS files and render merged report. Requires presence of the individual <sample_name>_QC.rds files present in the sample_sheet. [default: FALSE]\n",
    "  --skip_merge <TRUE/FALSE> Skip automatic merge after successful multi-sample render [default: FALSE]\n",
    "  --run_integration <TRUE/FALSE> Run integration and render integration report after successful merge [default: FALSE]\n",
    "  --integration_template <path> Path to integrate_RNA.Rmd [default: same dir as script]\n",
    "  --integration_only <TRUE/FALSE> Run only integration from existing merged_QC.rds [default: FALSE]\n",
    "  --integration_level <Batch|Sample> Required when integration runs\n",
    "  --help                  Show this help message\n",
    sep = ""
  )
  quit(status = status)
}

parse_cli_args <- function(args) {
  defaults <- list(
    sample_sheet = NULL,
    output_dir = "QC",
    template = NULL,
    merge_template = NULL,
    integration_template = NULL,
    integration_level = NULL,
    outs_subdir = "outs",
    ncores = 1L,
    sample_name = NULL,
    merge_only = "FALSE",
    skip_merge = "FALSE",
    run_integration = "FALSE",
    integration_only = "FALSE"
  )

  if (length(args) == 0L) return(defaults)

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (identical(key, "--help") || identical(key, "-h")) {
      print_usage_and_exit(0L)
    }
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }
    if (i == length(args)) {
      stop("Missing value for argument: ", key, call. = FALSE)
    }
    val <- args[[i + 1L]]
    if (startsWith(val, "--")) {
      stop("Missing value for argument: ", key, call. = FALSE)
    }

    name <- sub("^--", "", key)
    if (!name %in% names(defaults)) {
      stop("Unknown argument: ", key, call. = FALSE)
    }
    defaults[[name]] <- val
    i <- i + 2L
  }

  defaults$ncores <- suppressWarnings(as.integer(defaults$ncores))
  if (is.na(defaults$ncores) || defaults$ncores < 1L) {
    stop("--ncores must be a positive integer.", call. = FALSE)
  }

  parse_logical_flag <- function(x, flag_name) {
    x_chr <- tolower(trimws(as.character(x)))
    if (x_chr %in% c("true", "t", "1", "yes", "y")) return(TRUE)
    if (x_chr %in% c("false", "f", "0", "no", "n")) return(FALSE)
    stop(flag_name, " must be TRUE or FALSE.", call. = FALSE)
  }
  defaults$merge_only <- parse_logical_flag(defaults$merge_only, "--merge_only")
  defaults$skip_merge <- parse_logical_flag(defaults$skip_merge, "--skip_merge")
  defaults$run_integration <- parse_logical_flag(defaults$run_integration, "--run_integration")
  defaults$integration_only <- parse_logical_flag(defaults$integration_only, "--integration_only")

  if (isTRUE(defaults$integration_only)) {
    defaults$run_integration <- TRUE
  }

  normalize_integration_level <- function(x) {
    x_chr <- trimws(as.character(x))
    if (!nzchar(x_chr)) return(NA_character_)
    x_low <- tolower(x_chr)
    if (x_low == "batch") return("Batch")
    if (x_low == "sample") return("Sample")
    NA_character_
  }
  defaults$integration_level <- normalize_integration_level(defaults$integration_level)

  if (isTRUE(defaults$run_integration) && is.na(defaults$integration_level)) {
    stop("--integration_level is required when integration runs; use Batch or Sample.", call. = FALSE)
  }

  defaults
}

opt <- parse_cli_args(commandArgs(trailingOnly = TRUE))

script_arg <- commandArgs(trailingOnly = FALSE)[
  startsWith(commandArgs(trailingOnly = FALSE), "--file=")
][1L]
script_dir <- tryCatch(
  dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = FALSE)),
  error = function(e) getwd()
)

# ── Validate required arguments ────────────────────────────────────────────────
if (is.null(opt$sample_sheet)) {
  stop("--sample_sheet is required. See script header for usage.", call. = FALSE)
}

# Resolve/validate template path only when rendering is needed.
if (!isTRUE(opt$merge_only)) {
  # Resolve template path (default: same directory as this script)
  if (is.null(opt$template)) {
    opt$template <- file.path(script_dir, "sample_QC.Rmd")
  }

  if (!file.exists(opt$template)) {
    stop("Template not found: ", opt$template, call. = FALSE)
  }
}

if (is.null(opt$merge_template)) {
  opt$merge_template <- file.path(script_dir, "merge_analysis.Rmd")
}

if (is.null(opt$integration_template)) {
  opt$integration_template <- file.path(script_dir, "integrate_RNA.Rmd")
}

if (isTRUE(opt$run_integration) && !file.exists(opt$integration_template)) {
  stop("Integration template not found: ", opt$integration_template, call. = FALSE)
}

if (isTRUE(opt$run_integration) && isTRUE(opt$skip_merge) && !isTRUE(opt$integration_only)) {
  stop("--run_integration TRUE requires merge to run; do not combine with --skip_merge TRUE.", call. = FALSE)
}

if (isTRUE(opt$integration_only) && isTRUE(opt$merge_only)) {
  stop("--integration_only cannot be combined with --merge_only.", call. = FALSE)
}

if (isTRUE(opt$integration_only) && !is.null(opt$sample_name)) {
  stop("--integration_only cannot be combined with --sample_name.", call. = FALSE)
}

# ── Read and validate sample sheet ────────────────────────────────────────────
samples <- read.csv(opt$sample_sheet, stringsAsFactors = FALSE, strip.white = TRUE)

required_cols <- c("sample_name", "cellranger_dir")
missing_cols  <- setdiff(required_cols, colnames(samples))
if (length(missing_cols) > 0) {
  stop("Sample sheet is missing required column(s): ",
       paste(missing_cols, collapse = ", "), call. = FALSE)
}

if (anyDuplicated(samples$sample_name)) {
  stop("Duplicate sample_name entries found in sample sheet.", call. = FALSE)
}

# Optional batch column; default to "A" when missing or blank.
if (!"batch" %in% colnames(samples)) {
  samples$batch <- "A"
} else {
  samples$batch <- trimws(as.character(samples$batch))
  samples$batch[is.na(samples$batch) | samples$batch == ""] <- "A"
}

# ── Merge helper ──────────────────────────────────────────────────────────────
merge_sample_rds <- function(samples_df, output_dir) {
  if (nrow(samples_df) < 2L) {
    message("merge_only requested, but fewer than 2 samples were provided; skipping merge.")
    return(invisible(FALSE))
  }

  rds_paths <- file.path(output_dir, paste0(samples_df$sample_name, "_QC.rds"))
  found <- file.exists(rds_paths)
  if (!all(found)) {
    stop(
      "Cannot merge: missing RDS file(s) for sample(s): ",
      paste(samples_df$sample_name[!found], collapse = ", "),
      call. = FALSE
    )
  }

  suppressPackageStartupMessages(library(Seurat))
  message("Merging ", length(rds_paths), " existing sample RDS files...")

  obj_list <- lapply(seq_along(rds_paths), function(i) {
    message("  Loading: ", basename(rds_paths[[i]]))
    obj <- readRDS(rds_paths[[i]])
    obj$sample_name <- samples_df$sample_name[[i]]
    obj$batch <- samples_df$batch[[i]]
    obj
  })

  merged_obj <- merge(
    obj_list[[1L]],
    y = obj_list[-1L],
    add.cell.ids = samples_df$sample_name,
    project = "merged_QC",
    merge.data = TRUE
  )

  merged_rds <- file.path(output_dir, "merged_QC.rds")
  message("Saving merged object: ", merged_rds)
  message("  Cells: ", ncol(merged_obj), "  Features: ", nrow(merged_obj))
  saveRDS(merged_obj, file = merged_rds)
  message("Merge complete.")
  invisible(TRUE)
}

render_merge_rmd <- function(output_dir, merge_template) {
  if (!file.exists(merge_template)) {
    warning("Merge template not found; skipping merge report render: ", merge_template, call. = FALSE)
    return(invisible(FALSE))
  }

  merged_qc_rds <- file.path(output_dir, "merged_QC.rds")
  if (!file.exists(merged_qc_rds)) {
    warning("Merged object not found; skipping merge report render: ", merged_qc_rds, call. = FALSE)
    return(invisible(FALSE))
  }

  message("Rendering merge report from template: ", merge_template)
  out_html <- file.path(output_dir, "merge_analysis.html")

  ok <- tryCatch({
    rmarkdown::render(
      input          = merge_template,
      output_file    = out_html,
      knit_root_dir  = output_dir,
      envir          = new.env(parent = globalenv()),
      quiet          = TRUE
    )
    TRUE
  }, error = function(e) {
    warning("merge_analysis.Rmd render failed: ", conditionMessage(e), call. = FALSE)
    FALSE
  })

  if (isTRUE(ok)) {
    message("Merge report rendered: ", out_html)
  }
  invisible(ok)
}

render_integration_rmd <- function(output_dir, integration_template) {
  if (!file.exists(integration_template)) {
    warning("Integration template not found; skipping integration report render: ", integration_template, call. = FALSE)
    return(invisible(FALSE))
  }

  merged_qc_rds <- file.path(output_dir, "merged_QC.rds")
  if (!file.exists(merged_qc_rds)) {
    warning("Merged object not found; skipping integration report render: ", merged_qc_rds, call. = FALSE)
    return(invisible(FALSE))
  }

  message("Rendering integration report from template: ", integration_template)
  out_html <- file.path(output_dir, "integrate_RNA.html")

  ok <- tryCatch({
    rmarkdown::render(
      input          = integration_template,
      output_file    = out_html,
      params         = list(
        merged_rds_path    = merged_qc_rds,
        integrated_rds_path = file.path(output_dir, "integrated.rds"),
        output_dir         = output_dir,
        integration_level  = opt$integration_level
      ),
      knit_root_dir  = output_dir,
      envir          = new.env(parent = globalenv()),
      quiet          = TRUE
    )
    TRUE
  }, error = function(e) {
    warning("integrate_RNA.Rmd render failed: ", conditionMessage(e), call. = FALSE)
    FALSE
  })

  if (isTRUE(ok)) {
    message("Integration report rendered: ", out_html)
  }
  invisible(ok)
}

# ── integration_only mode: render integration report from existing merged_QC.rds
if (isTRUE(opt$integration_only)) {
  dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(opt$output_dir, "chunk_logs"), showWarnings = FALSE, recursive = TRUE)
  message(strrep("-", 60))
  render_integration_rmd(opt$output_dir, opt$integration_template)
  quit(status = 0)
}

# ── merge_only mode: do not resolve matrices or render; just merge existing RDS
if (isTRUE(opt$merge_only)) {
  if (!is.null(opt$sample_name)) {
    stop("--merge_only cannot be combined with --sample_name.", call. = FALSE)
  }
  if (isTRUE(opt$skip_merge)) {
    stop("--merge_only cannot be combined with --skip_merge TRUE.", call. = FALSE)
  }
  dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(opt$output_dir, "chunk_logs"), showWarnings = FALSE, recursive = TRUE)
  merge_ok <- merge_sample_rds(samples, opt$output_dir)
  if (isTRUE(merge_ok) && nrow(samples) > 1L) {
    message(strrep("-", 60))
    render_merge_rmd(opt$output_dir, opt$merge_template)
    if (isTRUE(opt$run_integration)) {
      message(strrep("-", 60))
      render_integration_rmd(opt$output_dir, opt$integration_template)
    }
  }
  quit(status = 0)
}

# ── Filter to a single sample if --sample_name was supplied ───────────────────
if (!is.null(opt$sample_name)) {
  if (isTRUE(opt$run_integration)) {
    stop("--run_integration TRUE cannot be combined with --sample_name because merge is not run.", call. = FALSE)
  }
  samples <- samples[samples$sample_name == opt$sample_name, ]
  if (nrow(samples) == 0) {
    stop("--sample_name '", opt$sample_name, "' not found in sample sheet.", call. = FALSE)
  }
}

# ── Optional threshold columns (NA means "use MAD default") ───────────────────
threshold_cols <- c("min_nCount_RNA", "max_nCount_RNA",
                    "min_nFeature_RNA", "max_nFeature_RNA",
                    "max_percent_mt", "min_malat1")

# Add any missing threshold columns as NA so downstream code is always uniform
for (col in threshold_cols) {
  if (!col %in% colnames(samples)) samples[[col]] <- NA_real_
}

# ── Resolve filtered/raw matrix directories with auto-detection ───────────────
is_10x_mex_dir <- function(path) {
  if (!dir.exists(path)) return(FALSE)
  files <- list.files(path)
  has_matrix <- any(grepl("^matrix\\.mtx(\\.gz)?$", files, ignore.case = TRUE))
  has_barcodes <- any(grepl("^barcodes\\.tsv(\\.gz)?$", files, ignore.case = TRUE))
  has_features <- any(grepl("^(features|genes)\\.tsv(\\.gz)?$", files, ignore.case = TRUE))
  has_matrix && has_barcodes && has_features
}

resolve_matrix_dir <- function(base_dir, kind = c("filtered", "raw")) {
  kind <- match.arg(kind)

  if (!dir.exists(base_dir)) return(NA_character_)

  preferred <- if (kind == "filtered") {
    c("filtered_feature_bc_matrix", "filtered_gene_bc_matrices")
  } else {
    c("raw_feature_bc_matrix", "raw_gene_bc_matrices")
  }

  for (nm in preferred) {
    p <- file.path(base_dir, nm)
    if (is_10x_mex_dir(p)) return(normalizePath(p, mustWork = TRUE))
  }

  child_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
  if (length(child_dirs) == 0) return(NA_character_)

  candidate <- child_dirs[
    grepl(kind, basename(child_dirs), ignore.case = TRUE) &
      vapply(child_dirs, is_10x_mex_dir, logical(1))
  ]

  if (length(candidate) == 1L) {
    return(normalizePath(candidate, mustWork = TRUE))
  }

  if (length(candidate) > 1L) {
    ranked <- order(
      !grepl("feature_bc_matrix", basename(candidate), ignore.case = TRUE),
      nchar(basename(candidate)),
      basename(candidate)
    )
    warning(
      "Multiple ", kind, " matrix directories detected in ", base_dir,
      "; using: ", basename(candidate[ranked[1L]]),
      call. = FALSE
    )
    return(normalizePath(candidate[ranked[1L]], mustWork = TRUE))
  }

  # Fallback: recursive search for nested layouts (e.g., per_sample_outs/*/count/sample_*_feature_bc_matrix)
  nested_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = TRUE)
  if (length(nested_dirs) > 0) {
    nested_candidate <- nested_dirs[
      grepl(kind, basename(nested_dirs), ignore.case = TRUE) &
        vapply(nested_dirs, is_10x_mex_dir, logical(1))
    ]
    if (length(nested_candidate) >= 1L) {
      ranked <- order(
        !grepl("feature_bc_matrix", basename(nested_candidate), ignore.case = TRUE),
        !grepl("sample_", basename(nested_candidate), ignore.case = TRUE),
        nchar(nested_candidate),
        nested_candidate
      )
      if (length(nested_candidate) > 1L) {
        warning(
          "Multiple nested ", kind, " matrix directories detected in ", base_dir,
          "; using: ", nested_candidate[ranked[1L]],
          call. = FALSE
        )
      }
      return(normalizePath(nested_candidate[ranked[1L]], mustWork = TRUE))
    }
  }

  NA_character_
}

candidate_base_dirs <- function(cellranger_dir, outs_subdir) {
  primary <- if (nzchar(outs_subdir)) file.path(cellranger_dir, outs_subdir) else cellranger_dir

  cands <- c(
    primary,
    if (basename(primary) == "outs") dirname(primary) else character(0),
    file.path(cellranger_dir, "outs"),
    file.path(cellranger_dir, "count"),
    file.path(cellranger_dir, "count", "outs")
  )

  unique(normalizePath(cands, winslash = "/", mustWork = FALSE))
}

resolve_from_candidates <- function(cellranger_dir, outs_subdir, kind = c("filtered", "raw")) {
  kind <- match.arg(kind)
  cands <- candidate_base_dirs(cellranger_dir, outs_subdir)
  for (base in cands) {
    p <- resolve_matrix_dir(base, kind = kind)
    if (!is.na(p)) return(p)
  }
  NA_character_
}

# ── Build per-sample path objects ──────────────────────────────────────────────
samples$search_bases <- vapply(
  samples$cellranger_dir,
  function(d) paste(candidate_base_dirs(d, opt$outs_subdir), collapse = " | "),
  FUN.VALUE = character(1)
)

samples$filtered_path <- vapply(
  samples$cellranger_dir,
  function(d) resolve_from_candidates(d, opt$outs_subdir, kind = "filtered"),
  FUN.VALUE = character(1)
)
samples$raw_path <- vapply(
  samples$cellranger_dir,
  function(d) resolve_from_candidates(d, opt$outs_subdir, kind = "raw"),
  FUN.VALUE = character(1)
)

# Report resolved matrix directories for transparency/debugging in job logs
message("Resolved matrix directories:")
for (i in seq_len(nrow(samples))) {
  message("  ", samples$sample_name[[i]])
  message("    filtered: ", samples$filtered_path[[i]])
  message("    raw     : ", samples$raw_path[[i]])
}

# Warn if any input directories are missing (don't abort; let render report the error)
filtered_ok <- !is.na(samples$filtered_path) & dir.exists(samples$filtered_path)
raw_ok <- !is.na(samples$raw_path) & dir.exists(samples$raw_path)
missing_inputs <- samples[!(filtered_ok & raw_ok), c("sample_name", "search_bases")]
if (nrow(missing_inputs) > 0) {
  detail <- apply(missing_inputs, 1, function(x) {
    paste0(x[["sample_name"]], " (searched bases: ", x[["search_bases"]], ")")
  })
  warning(
    "Could not resolve filtered/raw matrix directories for:\n  ",
    paste(detail, collapse = "\n  "),
    call. = FALSE
  )
}

# ── Create output directories ─────────────────────────────────────────────────
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(opt$output_dir, "chunk_logs"), showWarnings = FALSE, recursive = TRUE)

# ── Render function for a single sample ───────────────────────────────────────
render_sample <- function(row) {
  sample_name <- row[["sample_name"]]
  message("[", Sys.time(), "] Starting: ", sample_name)
  message("[", Sys.time(), "]   filtered_path: ", row[["filtered_path"]])
  message("[", Sys.time(), "]   raw_path     : ", row[["raw_path"]])

  out_html <- file.path(opt$output_dir,
                        paste0(sample_name, "_QC.html"))

  result <- tryCatch({
    rmarkdown::render(
      input         = opt$template,
      output_file   = out_html,
      params        = list(
        sample_name      = sample_name,
        batch            = row[["batch"]],
        filtered_path    = row[["filtered_path"]],
        raw_path         = row[["raw_path"]],
        output_dir       = opt$output_dir,
        min_nCount_RNA   = suppressWarnings(as.numeric(row[["min_nCount_RNA"]])),
        max_nCount_RNA   = suppressWarnings(as.numeric(row[["max_nCount_RNA"]])),
        min_nFeature_RNA = suppressWarnings(as.numeric(row[["min_nFeature_RNA"]])),
        max_nFeature_RNA = suppressWarnings(as.numeric(row[["max_nFeature_RNA"]])),
        max_percent_mt   = suppressWarnings(as.numeric(row[["max_percent_mt"]])),
        min_malat1       = suppressWarnings(as.numeric(row[["min_malat1"]]))
      ),
      envir         = new.env(parent = globalenv()),
      quiet         = TRUE
    )
    "SUCCESS"
  }, error = function(e) {
    paste("FAILED:", conditionMessage(e))
  })

  message("[", Sys.time(), "] ", sample_name, " — ", result)
  data.frame(sample_name = sample_name, status = result,
             stringsAsFactors = FALSE)
}

# ── Render all samples (sequential or parallel) ────────────────────────────────
message("Rendering ", nrow(samples), " sample(s) using ", opt$ncores, " core(s)...")
message("Template:   ", opt$template)
message("Output dir: ", opt$output_dir)
message(strrep("-", 60))

sample_list <- split(samples, seq_len(nrow(samples)))

if (opt$ncores > 1L) {
  # mclapply forks child processes; each gets its own R session
  results_list <- parallel::mclapply(sample_list, render_sample,
                                     mc.cores = opt$ncores)
} else {
  results_list <- lapply(sample_list, render_sample)
}

# ── Summary report ─────────────────────────────────────────────────────────────
results <- do.call(rbind, results_list)
message(strrep("-", 60))
message("SUMMARY:")
print(results, row.names = FALSE)

failed <- results[!grepl("^SUCCESS", results$status), ]
if (nrow(failed) > 0) {
  message("\nFailed samples:")
  print(failed, row.names = FALSE)
  quit(status = 1)
} else {
  message("\nAll samples completed successfully.")
}

# ── Merge all per-sample Seurat objects when more than one sample was rendered ─
# Only attempt merge when processing the full sample sheet (not a single-sample run)
# and only when every sample succeeded.
if (is.null(opt$sample_name) && nrow(samples) > 1L && nrow(failed) == 0L) {
  if (isTRUE(opt$skip_merge)) {
    message(strrep("-", 60))
    message("Skipping merge and merge report because --skip_merge is TRUE.")
  } else {
    message(strrep("-", 60))
    merge_ok <- merge_sample_rds(samples, opt$output_dir)
    if (isTRUE(merge_ok)) {
      message(strrep("-", 60))
      render_merge_rmd(opt$output_dir, opt$merge_template)
      if (isTRUE(opt$run_integration)) {
        message(strrep("-", 60))
        render_integration_rmd(opt$output_dir, opt$integration_template)
      }
    }
  }
}
