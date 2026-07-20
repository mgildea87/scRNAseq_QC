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
#   sample_name,filtered_path,raw_path,batch
#   Control_1,/path/to/filtered_feature_bc_matrix,/path/to/raw_feature_bc_matrix,A
#   Treatment_1,/path/to/filtered_feature_bc_matrix,/path/to/raw_feature_bc_matrix,A
#
# The pipeline now expects explicit matrix paths in the sample sheet.
# Each path should point directly to a 10x MEX directory or a 10x .h5/.hdf5 file.
#
# OPTIONAL ARGUMENTS:
#   --template      Path to the .Rmd template
#                   (default: same directory as this script)
#   --merge_template Path to merge_analysis.Rmd template
#                   (default: same directory as this script)
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
#   --use_cellbender If TRUE, prefer/require cellbender_filtered.h5 for
#                   filtered input resolution (default: FALSE)
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
    "  --ncores <int>          Number of samples to render in parallel [default: 1]\n",
    "  --sample_name <name>    Process only one sample_name from sheet\n",
    "  --merge_only <TRUE/FALSE>  If TRUE, skip individual sample QC and only merge existing RDS files and render merged report. Requires presence of the individual <sample_name>_QC.rds files present in the sample_sheet. [default: FALSE]\n",
    "  --skip_merge <TRUE/FALSE> Skip automatic merge after successful multi-sample render [default: FALSE]\n",
    "  --run_integration <TRUE/FALSE> Run integration and render integration report after successful merge [default: FALSE]\n",
    "  --integration_template <path> Path to integrate_RNA.Rmd [default: same dir as script]\n",
    "  --integration_only <TRUE/FALSE> Run only integration from existing merged_QC.rds [default: FALSE]\n",
    "  --integration_level <Batch|Sample> Required when integration runs\n",
    "  --use_cellbender <TRUE/FALSE> Use cellbender_filtered.h5 as filtered input [default: FALSE]\n",
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
    ncores = 1L,
    sample_name = NULL,
    merge_only = "FALSE",
    skip_merge = "FALSE",
    run_integration = "FALSE",
    integration_only = "FALSE",
    use_cellbender = "FALSE"
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
  defaults$use_cellbender <- parse_logical_flag(defaults$use_cellbender, "--use_cellbender")

  if (isTRUE(defaults$integration_only)) {
    defaults$run_integration <- TRUE
  }

  normalize_integration_level <- function(x) {
    if (is.null(x) || length(x) == 0L) return(NA_character_)
    x_chr <- trimws(as.character(x[[1L]]))
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

message("Run options:")
message("  sample_sheet    : ", opt$sample_sheet)
message("  output_dir      : ", opt$output_dir)
message("  sample_name     : ", if (is.null(opt$sample_name)) "<all>" else opt$sample_name)
message("  ncores          : ", opt$ncores)
message("  use_cellbender  : ", opt$use_cellbender)
message("  skip_merge      : ", opt$skip_merge)
message("  run_integration : ", opt$run_integration)
message("  integration_only: ", opt$integration_only)
message("  integration_level: ", if (is.na(opt$integration_level)) "<unset>" else opt$integration_level)

script_arg <- commandArgs(trailingOnly = FALSE)[
  startsWith(commandArgs(trailingOnly = FALSE), "--file=")
][1L]
script_dir <- tryCatch(
  dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = FALSE)),
  error = function(e) getwd()
)

# Normalize output path early so metadata and outputs are colocated.
opt$output_dir <- normalizePath(opt$output_dir, winslash = "/", mustWork = FALSE)

safe_system_output <- function(cmd, args = character()) {
  out <- tryCatch(
    suppressWarnings(system2(cmd, args = args, stdout = TRUE, stderr = TRUE)),
    error = function(e) character()
  )
  if (length(out) == 0L) NA_character_ else paste(out, collapse = "\n")
}

write_run_metadata <- function(opt, script_dir) {
  dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
  meta_dir <- file.path(opt$output_dir, "run_metadata")
  dir.create(meta_dir, showWarnings = FALSE, recursive = TRUE)

  invocation_id <- paste0(format(Sys.time(), "%Y%m%dT%H%M%S"), "_", Sys.getpid())

  git_root <- safe_system_output("git", c("-C", script_dir, "rev-parse", "--show-toplevel"))
  if (is.na(git_root)) git_root <- script_dir

  repo_url <- safe_system_output("git", c("-C", git_root, "remote", "get-url", "origin"))
  commit_sha <- safe_system_output("git", c("-C", git_root, "rev-parse", "HEAD"))
  branch_name <- safe_system_output("git", c("-C", git_root, "rev-parse", "--abbrev-ref", "HEAD"))
  exact_tag <- safe_system_output("git", c("-C", git_root, "describe", "--tags", "--exact-match"))
  describe_tag <- safe_system_output("git", c("-C", git_root, "describe", "--tags", "--always"))
  worktree_status <- safe_system_output("git", c("-C", git_root, "status", "--porcelain"))
  worktree_label <- if (is.na(worktree_status) || worktree_status == "") "clean" else "dirty"

  release_path <- file.path(meta_dir, "release_info.txt")
  writeLines(c(
    paste0("recorded_utc=", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("git_root=", git_root),
    paste0("repo_url=", if (is.na(repo_url)) "<unknown>" else repo_url),
    paste0("commit_sha=", if (is.na(commit_sha)) "<unknown>" else commit_sha),
    paste0("branch=", if (is.na(branch_name)) "<unknown>" else branch_name),
    paste0("exact_tag=", if (is.na(exact_tag)) "<none>" else exact_tag),
    paste0("describe=", if (is.na(describe_tag)) "<none>" else describe_tag),
    paste0("worktree_status=", worktree_label)
  ), con = release_path)

  trailing <- commandArgs(trailingOnly = TRUE)
  run_cmd <- paste(c("Rscript", shQuote(file.path(script_dir, "qc_batch_runner.R")), shQuote(trailing)), collapse = " ")
  run_cmd_path <- file.path(meta_dir, paste0("run_command_", invocation_id, ".sh"))
  writeLines(c("#!/bin/bash", "set -euo pipefail", run_cmd), con = run_cmd_path)

  conda_info_path <- file.path(meta_dir, "conda_info.txt")
  conda_env <- Sys.getenv("CONDA_DEFAULT_ENV", unset = "")
  conda_envs <- safe_system_output("conda", c("info", "--envs"))
  writeLines(c(
    paste0("recorded_utc=", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("CONDA_DEFAULT_ENV=", ifelse(nzchar(conda_env), conda_env, "<unset>")),
    "",
    "conda info --envs:",
    if (is.na(conda_envs)) "<conda unavailable>" else conda_envs
  ), con = conda_info_path)

  if (nzchar(conda_env) && !is.na(conda_envs)) {
    conda_export <- safe_system_output("conda", c("env", "export", "-n", conda_env))
    if (!is.na(conda_export)) {
      writeLines(conda_export, con = file.path(meta_dir, "conda_env_export.yml"))
    }
  }

  start_path <- file.path(meta_dir, paste0("run_start_", invocation_id, ".txt"))
  writeLines(c(
    paste0("start_utc=", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("invocation_id=", invocation_id),
    paste0("command=", run_cmd)
  ), con = start_path)

  list(meta_dir = meta_dir, invocation_id = invocation_id, run_cmd = run_cmd)
}

run_meta <- write_run_metadata(opt, script_dir)
run_status <- "FAILED"
run_exit_status <- 1L

on.exit({
  end_path <- file.path(run_meta$meta_dir, paste0("run_end_", run_meta$invocation_id, ".txt"))
  writeLines(c(
    paste0("end_utc=", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("invocation_id=", run_meta$invocation_id),
    paste0("status=", run_status),
    paste0("exit_status=", run_exit_status),
    paste0("command=", run_meta$run_cmd)
  ), con = end_path)
}, add = TRUE)

finish_run <- function(status_code = 0L, status_label = "SUCCESS") {
  run_status <<- status_label
  run_exit_status <<- as.integer(status_code)
}

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

required_cols <- c("sample_name", "filtered_path", "raw_path")
missing_cols  <- setdiff(required_cols, colnames(samples))
if (length(missing_cols) > 0) {
  stop(
    "Sample sheet is missing required column(s): ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
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

normalize_optional_path <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == ""] <- NA_character_
  x
}

samples$filtered_path <- normalize_optional_path(samples$filtered_path)
samples$raw_path <- normalize_optional_path(samples$raw_path)

if (anyNA(samples$filtered_path) || anyNA(samples$raw_path)) {
  stop(
    "Every sample must provide non-empty filtered_path and raw_path values.",
    call. = FALSE
  )
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
  finish_run(0L, "SUCCESS")
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
  finish_run(0L, "SUCCESS")
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

# ── Build per-sample path objects ──────────────────────────────────────────────
samples$filtered_path <- samples$filtered_path
samples$raw_path <- samples$raw_path

message("Matrix inputs:")
for (i in seq_len(nrow(samples))) {
  message("  ", samples$sample_name[[i]])
  message("    filtered: ", samples$filtered_path[[i]])
  message("    raw     : ", samples$raw_path[[i]])
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
    template_support_dir <- file.path(dirname(normalizePath(opt$template, mustWork = TRUE)), "support_files")
    rmarkdown::render(
      input         = opt$template,
      output_file   = out_html,
      params        = list(
        sample_name      = sample_name,
        batch            = row[["batch"]],
        filtered_path    = row[["filtered_path"]],
        raw_path         = row[["raw_path"]],
        use_cellbender   = isTRUE(opt$use_cellbender),
        support_dir      = template_support_dir,
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
  finish_run(1L, "FAILED")
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

finish_run(0L, "SUCCESS")
