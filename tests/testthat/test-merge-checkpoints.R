source(file.path(repo_root, "support_files", "merge_checkpoint_helpers.R"))

testthat::test_that("LISI checkpoints are reused only for matching inputs and results", {
  checkpoint_dir <- tempfile("lisi-checkpoint-")
  dir.create(checkpoint_dir)
  path <- file.path(checkpoint_dir, "lisi.rds")
  calls <- 0L
  compute <- function() {
    calls <<- calls + 1L
    list(RNA = data.frame(batch = calls), SCT = data.frame(batch = calls))
  }
  inputs <- list(pca = matrix(1:4, nrow = 2), metadata = c("A", "B"), perplexity = 30)
  fingerprint <- checkpoint_fingerprint(inputs)

  original <- get_lisi_checkpoint(path, fingerprint, FALSE, compute)
  testthat::expect_identical(readRDS(path), original)
  testthat::expect_identical(get_lisi_checkpoint(path, fingerprint, TRUE, compute), original)
  testthat::expect_equal(calls, 1L)

  inputs$metadata[[1L]] <- "C"
  changed <- checkpoint_fingerprint(inputs)
  get_lisi_checkpoint(path, changed, TRUE, compute)
  testthat::expect_equal(calls, 2L)
  inputs$perplexity <- 40
  changed <- checkpoint_fingerprint(inputs)
  get_lisi_checkpoint(path, changed, TRUE, compute)
  testthat::expect_equal(calls, 3L)

  unlink(paste0(path, ".provenance.rds"))
  get_lisi_checkpoint(path, changed, TRUE, compute)
  testthat::expect_equal(calls, 4L)

  saveRDS(list(incomplete = TRUE), path)
  get_lisi_checkpoint(path, changed, TRUE, compute)
  testthat::expect_equal(calls, 5L)

  get_lisi_checkpoint(path, changed, FALSE, compute)
  testthat::expect_equal(calls, 6L)

  saveRDS(list(version = 1L), paste0(path, ".provenance.rds"))
  testthat::expect_error(get_lisi_checkpoint(path, changed, TRUE, compute),
                         "Invalid LISI checkpoint provenance")

  writeLines("not an RDS", path)
  saveRDS(list(version = 1L, input_hash = changed, result_hash = checkpoint_file_hash(path)),
          paste0(path, ".provenance.rds"))
  testthat::expect_error(get_lisi_checkpoint(path, changed, TRUE, compute))
})

testthat::test_that("merged object is saved only when preprocessing changes it", {
  testthat::skip_if_not_installed("Seurat")
  set.seed(1)
  counts <- matrix(rpois(40 * 12, 2), nrow = 40,
                   dimnames = list(paste0("g", 1:40), paste0("c", 1:12)))
  object <- Seurat::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  object$batch <- rep(c("A", "B"), each = 6)
  object$sample_name <- object$batch
  object$library_id <- object$batch
  object$marker_sample_id <- object$batch
  object$RNA_merged_clusters <- rep("0", 12)
  object$SCT_merged_clusters <- rep("0", 12)
  for (reduction in c("PCA_RNA_merged", "UMAP_RNA_merged",
                      "PCA_SCT_merged", "UMAP_SCT_merged")) {
    key <- paste0(tolower(gsub("_", "", reduction)), "_")
    embeddings <- matrix(rnorm(12 * 2), nrow = 12,
                         dimnames = list(colnames(object), paste0(key, 1:2)))
    object[[reduction]] <- Seurat::CreateDimReducObject(
      embeddings = embeddings, key = key, assay = "RNA"
    )
  }
  object[["RNA"]] <- split(object[["RNA"]], f = object$batch)

  checkpoint_dir <- tempfile("merged-checkpoint-")
  dir.create(checkpoint_dir)
  withr::local_dir(checkpoint_dir)
  saveRDS(object, "merged_QC.rds")
  lines <- readLines(file.path(repo_root, "merge_analysis.Rmd"), warn = FALSE)
  start <- grep("^```\\{r load-merged-object", lines)
  end <- start + which(lines[(start + 1L):length(lines)] == "```")[[1L]]
  code <- paste(lines[(start + 1L):(end - 1L)], collapse = "\n")
  env <- new.env(parent = asNamespace("Seurat"))
  env$rna_layers_need_join <- rna_layers_need_join
  env$validate_marker_groups <- validate_marker_groups
  writes <- 0L
  env$atomic_save_rds <- function(object, path) {
    writes <<- writes + 1L
    atomic_save_rds(object, path)
  }

  old_object <- object
  old_object$marker_sample_id <- NULL
  saveRDS(old_object, "merged_QC.rds")
  testthat::expect_error(eval(parse(text = code), env), "Re-merge the per-sample QC objects")
  saveRDS(object, "merged_QC.rds")
  eval(parse(text = code), env)
  testthat::expect_equal(writes, 1L)
  testthat::expect_false(rna_layers_need_join(readRDS("merged_QC.rds")[["RNA"]]))
  hash <- checkpoint_file_hash("merged_QC.rds")
  eval(parse(text = code), env)
  testthat::expect_equal(writes, 1L)
  testthat::expect_identical(checkpoint_file_hash("merged_QC.rds"), hash)

  incomplete <- readRDS("merged_QC.rds")
  pca <- incomplete[["PCA_RNA_merged"]]
  incomplete[["PCA_RNA_merged"]] <- NULL
  saveRDS(incomplete, "merged_QC.rds")
  env$NormalizeData <- env$FindVariableFeatures <- env$ScaleData <-
    env$RunUMAP <- env$FindNeighbors <- env$FindClusters <- function(object, ...) object
  env$RunPCA <- function(object, ...) {
    object[["PCA_RNA_merged"]] <- pca
    object
  }
  eval(parse(text = code), env)
  testthat::expect_equal(writes, 2L)
  testthat::expect_true("PCA_RNA_merged" %in% names(readRDS("merged_QC.rds")@reductions))
})

testthat::test_that("marker grouping rejects incomplete or ambiguous identities", {
  metadata <- data.frame(HTO_sample = c("tag1", "tag2"),
                         replicate_id = c("A_tag1", "A_tag2"))
  testthat::expect_identical(marker_groups_for_library(metadata, "A", TRUE),
                             c("A_tag1", "A_tag2"))
  testthat::expect_identical(marker_groups_for_library(metadata, "B", FALSE), c("B", "B"))
  testthat::expect_error(marker_groups_for_library(metadata["HTO_sample"], "A", TRUE),
                         "missing metadata")
  metadata$replicate_id[[2L]] <- NA_character_
  testthat::expect_error(marker_groups_for_library(metadata, "A", TRUE), "Invalid HTO_sample or replicate_id")
  metadata$replicate_id[[2L]] <- "B_tag2"
  testthat::expect_error(marker_groups_for_library(metadata, "A", TRUE), "Invalid HTO_sample or replicate_id")
  testthat::expect_error(validate_marker_groups(data.frame(sample_name = "A")),
                         "Re-merge the per-sample QC objects")
  testthat::expect_error(validate_marker_groups(data.frame(marker_sample_id = c("tag", "tag"),
                                                           library_id = c("A", "B"))),
                         "shared across libraries")
  testthat::expect_error(validate_marker_groups(data.frame(marker_sample_id = NA_character_,
                                                           library_id = "A")), "missing or blank")
})

testthat::test_that("ordinary and merge-only runs retain distinct HTO replicates", {
  testthat::skip_if_not_installed("Seurat")
  testthat::skip_if_not_installed("rmarkdown")
  testthat::skip_if_not(rmarkdown::pandoc_available())
  work <- tempfile("marker-merge-")
  dir.create(work)
  sheet_path <- file.path(work, "samples.csv")
  write_transposed_sample_sheet(sheet_path, c("A", "B", "C"))
  sheet <- read.csv(sheet_path, colClasses = "character", check.names = FALSE)
  sheet[sheet$field == "use_hashtag", c("A", "B", "C")] <- c("TRUE", "TRUE", "FALSE")
  write.csv(sheet, sheet_path, row.names = FALSE, na = "NA")

  sample_template <- file.path(work, "sample.Rmd")
  param_names <- names(rmarkdown::yaml_front_matter(file.path(repo_root, "sample_QC.Rmd"))$params)
  writeLines(c("---", "output: html_document", "params:",
               paste0("  ", param_names, ": ''"), "---",
               "```{r}",
               "stopifnot(identical(params$use_hashtag, params$sample_name != 'C'))",
               "```"), sample_template)
  merge_template <- file.path(work, "merge.Rmd")
  writeLines(c("---", "output: html_document", "---", "Merged."), merge_template)
  runner <- file.path(repo_root, "qc_batch_runner.R")

  for (merge_only in c(TRUE, FALSE)) {
    output_dir <- file.path(work, if (merge_only) "merge-only" else "ordinary")
    dir.create(output_dir)
    for (library in c("A", "B", "C")) {
      counts <- Matrix::Matrix(matrix(1, nrow = 2, ncol = 2,
                                      dimnames = list(c("g1", "g2"), c("cell1", "cell2"))),
                               sparse = TRUE)
      object <- Seurat::CreateSeuratObject(counts)
      if (library != "C") {
        object$HTO_sample <- c("tag1", "tag2")
        object$replicate_id <- paste(library, object$HTO_sample, sep = "_")
      }
      saveRDS(object, file.path(output_dir, paste0(library, "_QC.rds")))
    }
    args <- c("--vanilla", shQuote(runner), "--sample_sheet", shQuote(sheet_path),
              "--output_dir", shQuote(output_dir), "--template", shQuote(sample_template),
              "--merge_template", shQuote(merge_template),
              "--merge_only", if (merge_only) "TRUE" else "FALSE")
    output <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), args,
                                       stdout = TRUE, stderr = TRUE))
    status <- attr(output, "status")
    testthat::expect_true(is.null(status) || status == 0L, info = paste(output, collapse = "\n"))
    testthat::expect_true(file.exists(file.path(output_dir, "merge_analysis.html")))
    merged <- readRDS(file.path(output_dir, "merged_QC.rds"))
    metadata <- merged@meta.data
    testthat::expect_identical(sort(unique(metadata$marker_sample_id)),
                               c("A_tag1", "A_tag2", "B_tag1", "B_tag2", "C"))
    testthat::expect_identical(sort(unique(metadata$sample_name[metadata$library_id == "A"])),
                               c("tag1", "tag2"))
    testthat::expect_true(all(metadata$sample_name[metadata$library_id == "C"] == "C"))
  }

  sheet[sheet$field == "use_hashtag", "B"] <- "invalid"
  write.csv(sheet, sheet_path, row.names = FALSE, na = "NA")
  output <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
                                     c("--vanilla", shQuote(runner), "--sample_sheet", shQuote(sheet_path),
                                       "--output_dir", shQuote(file.path(work, "invalid")),
                                       "--merge_only", "TRUE"),
                                     stdout = TRUE, stderr = TRUE))
  testthat::expect_false(is.null(attr(output, "status")))
  testthat::expect_match(paste(output, collapse = "\n"), "Invalid use_hashtag value")
})

testthat::test_that("all marker reports use the marker sample grouping", {
  for (file in c("merge_analysis.Rmd", "integrate_RNA.Rmd")) {
    lines <- readLines(file.path(repo_root, file), warn = FALSE)
    calls <- grep("^FindMarkersBulk\\(", lines, value = TRUE)
    testthat::expect_length(calls, if (file == "merge_analysis.Rmd") 2L else 4L)
    testthat::expect_true(all(grepl("sample_ident = 'marker_sample_id'", calls, fixed = TRUE)))
  }
})

testthat::test_that("legacy checkpoints fail both report resume paths", {
  for (package in c("Seurat", "rmarkdown", "lisi", "CVRCFunc", "harmony", "pheatmap")) {
    testthat::skip_if_not_installed(package)
  }
  testthat::skip_if_not(rmarkdown::pandoc_available())
  work <- tempfile("legacy-marker-checkpoint-")
  dir.create(work)
  dir.create(file.path(work, "support_files"))
  for (file in c("qc_helpers.R", "merge_checkpoint_helpers.R")) {
    file.copy(file.path(repo_root, "support_files", file), file.path(work, "support_files", file))
  }
  for (file in c("merge_analysis.Rmd", "integrate_RNA.Rmd")) {
    file.copy(file.path(repo_root, file), file.path(work, file))
  }
  sheet <- file.path(work, "samples.csv")
  write_transposed_sample_sheet(sheet)
  counts <- Matrix::Matrix(matrix(1, nrow = 2, ncol = 2,
                                  dimnames = list(c("g1", "g2"), c("cell1", "cell2"))),
                           sparse = TRUE)
  object <- Seurat::CreateSeuratObject(counts)
  object$sample_name <- "A"
  object$library_id <- "A"
  object$batch <- "A"
  saveRDS(object, file.path(work, "merged_QC.rds"))
  runner <- file.path(repo_root, "qc_batch_runner.R")
  common <- c("--vanilla", shQuote(runner), "--sample_sheet", shQuote(sheet),
              "--output_dir", shQuote(work))
  for (args in list(c("--resume_merge_report", "TRUE", "--merge_template",
                       shQuote(file.path(work, "merge_analysis.Rmd"))),
                    c("--integration_only", "TRUE", "--integration_template",
                      shQuote(file.path(work, "integrate_RNA.Rmd"))))) {
    output <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), c(common, args),
                                       stdout = TRUE, stderr = TRUE))
    testthat::expect_false(is.null(attr(output, "status")), info = paste(output, collapse = "\n"))
    testthat::expect_match(paste(output, collapse = "\n"), "Re-merge the per-sample QC objects")
  }
})

testthat::test_that("the report computes and resumes real LISI results", {
  testthat::skip_if_not_installed("Seurat")
  testthat::skip_if_not_installed("lisi")
  set.seed(2)
  counts <- matrix(rpois(40 * 100, 2), nrow = 40,
                   dimnames = list(paste0("g", 1:40), paste0("c", 1:100)))
  object <- Seurat::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  object$batch <- rep(c("A", "B"), each = 50)
  object$sample_name <- rep(c("S1", "S2", "S3", "S4"), each = 25)
  for (reduction in c("PCA_RNA_merged", "PCA_SCT_merged")) {
    key <- paste0(tolower(gsub("_", "", reduction)), "_")
    embeddings <- matrix(rnorm(100 * 30), nrow = 100,
                         dimnames = list(colnames(object), paste0(key, 1:30)))
    object[[reduction]] <- Seurat::CreateDimReducObject(
      embeddings = embeddings, key = key, assay = "RNA"
    )
  }
  checkpoint_dir <- tempfile("real-lisi-")
  dir.create(checkpoint_dir)
  withr::local_dir(checkpoint_dir)
  lines <- readLines(file.path(repo_root, "merge_analysis.Rmd"), warn = FALSE)
  start <- grep("^```\\{r compute-lisi", lines)
  end <- start + which(lines[(start + 1L):length(lines)] == "```")[[1L]]
  code <- paste(lines[(start + 1L):(end - 1L)], collapse = "\n")
  env <- new.env(parent = globalenv())
  env$seurat <- object
  env$params <- list(resume_lisi = FALSE)
  env$checkpoint_fingerprint <- checkpoint_fingerprint
  env$get_lisi_checkpoint <- get_lisi_checkpoint
  calls <- 0L
  env$compute_lisi <- function(...) {
    calls <<- calls + 1L
    lisi::compute_lisi(...)
  }
  eval(parse(text = code), env)
  original <- env$lisi
  testthat::expect_equal(calls, 2L)
  testthat::expect_identical(names(original), c("RNA", "SCT"))
  testthat::expect_identical(dim(original$RNA), c(100L, 2L))
  testthat::expect_identical(colnames(original$SCT), c("batch", "sample_name"))
  env$params$resume_lisi <- TRUE
  eval(parse(text = code), env)
  testthat::expect_equal(calls, 2L)
  testthat::expect_identical(env$lisi, original)
  env$seurat$batch[[1L]] <- "B"
  eval(parse(text = code), env)
  testthat::expect_equal(calls, 4L)
})

testthat::test_that("resume runner renders without merging and reports failures", {
  testthat::skip_if_not_installed("rmarkdown")
  testthat::skip_if_not(rmarkdown::pandoc_available())
  work <- tempfile("resume-runner-")
  dir.create(work)
  sample_sheet <- file.path(work, "samples.csv")
  write_transposed_sample_sheet(sample_sheet)
  merged <- file.path(work, "merged_QC.rds")
  saveRDS(list(unchanged = TRUE), merged)
  original_hash <- checkpoint_file_hash(merged)
  template <- file.path(work, "report.Rmd")
  writeLines(c(
    "---", "output: html_document", "params:", "  resume_lisi: false", "---",
    "```{r}", "writeLines(as.character(params$resume_lisi), 'resume-observed.txt')", "```"
  ), template)
  runner <- file.path(repo_root, "qc_batch_runner.R")
  args <- c(shQuote(runner), "--sample_sheet", shQuote(sample_sheet),
            "--output_dir", shQuote(work), "--merge_template", shQuote(template),
            "--resume_merge_report", "TRUE")
  run <- function(args) {
    output <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
                                      c("--vanilla", args), stdout = TRUE, stderr = TRUE))
    status <- attr(output, "status")
    list(output = output, status = if (is.null(status)) 0L else status)
  }
  result <- run(args)
  testthat::expect_equal(result$status, 0L, info = paste(result$output, collapse = "\n"))
  testthat::expect_identical(readLines(file.path(work, "resume-observed.txt")), "TRUE")
  testthat::expect_identical(checkpoint_file_hash(merged), original_hash)
  testthat::expect_true(file.exists(file.path(work, "merge_analysis.html")))

  writeLines(c(
    "---", "output: html_document", "params:", "  resume_lisi: false", "---",
    "```{r}", "stop('intentional render failure')", "```"
  ), template)
  result <- run(args)
  testthat::expect_equal(result$status, 1L)
  testthat::expect_match(paste(result$output, collapse = "\n"), "intentional render failure")
  result <- run(c(args, "--merge_only", "TRUE"))
  testthat::expect_false(result$status == 0L)
  testthat::expect_match(paste(result$output, collapse = "\n"), "cannot be combined")
  testthat::expect_identical(checkpoint_file_hash(merged), original_hash)
  unlink(merged)
  result <- run(args)
  testthat::expect_equal(result$status, 1L)
  testthat::expect_match(paste(result$output, collapse = "\n"), "Merged object not found")
})
