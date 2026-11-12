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
  writes <- 0L
  env$atomic_save_rds <- function(object, path) {
    writes <<- writes + 1L
    atomic_save_rds(object, path)
  }

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
