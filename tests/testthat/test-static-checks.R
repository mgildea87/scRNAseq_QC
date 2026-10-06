testthat::test_that("R scripts and R Markdown parse", {
  testthat::expect_no_error(invisible(parse(file.path(repo_root, "qc_batch_runner.R"))))

  generated_r <- tempfile(fileext = ".R")
  on.exit(unlink(generated_r), add = TRUE)
  template_path <- file.path(repo_root, "sample_QC.Rmd")
  purl_env <- new.env(parent = globalenv())
  purl_env$params <- rmarkdown::yaml_front_matter(template_path)$params
  knitr::purl(template_path, output = generated_r, quiet = TRUE, envir = purl_env)
  testthat::expect_no_error(invisible(parse(generated_r)))

  template_text <- paste(readLines(file.path(repo_root, "sample_QC.Rmd"), warn = FALSE), collapse = "\n")
  testthat::expect_match(template_text, "positive.quantile = 0.99", fixed = TRUE)
  testthat::expect_match(template_text, "hto_singlet_keep_mask", fixed = TRUE)
  testthat::expect_match(template_text, "hto_data <- GetAssayData(seurat, assay = 'HTO', layer = 'data')", fixed = TRUE)
  testthat::expect_match(template_text, "hto_singlet_expression_thresholds", fixed = TRUE)
  testthat::expect_match(template_text, "requireNamespace('scDblFinder'", fixed = TRUE)
  testthat::expect_match(template_text, "scDblFinder::recoverDoublets", fixed = TRUE)
  testthat::expect_match(template_text, "reducedDimNames(hto_sce)", fixed = TRUE)
  testthat::expect_match(template_text, "use.dimred = hto_dimreds$pca", fixed = TRUE)
  testthat::expect_match(template_text, "dimred = hto_dimreds$umap", fixed = TRUE)
  testthat::expect_no_match(template_text, "scran::recoverDoublets", fixed = TRUE)
  testthat::expect_no_match(template_text, "requireNamespace('scran'", fixed = TRUE)
  testthat::expect_match(template_text, "validate_raw_adt_path(params$use_adt, params$raw_adt_path)", fixed = TRUE)
  testthat::expect_match(template_text, "raw_adt_path <- ''", fixed = TRUE)
  testthat::expect_match(template_text, "input_path = params$raw_path, label = 'raw'", fixed = TRUE)
  testthat::expect_match(template_text, "adt_path = raw_adt_path", fixed = TRUE)
  testthat::expect_match(template_text, "geom_vline(xintercept = hto_threshold, color = 'red'", fixed = TRUE)
  testthat::expect_no_match(template_text, "parse_hto_demux_thresholds", fixed = TRUE)
  testthat::expect_no_match(template_text, "params$hto_positive_quantile", fixed = TRUE)
  testthat::expect_no_match(template_text, "params$hto_filter_doublets_negatives", fixed = TRUE)
})

testthat::test_that("recoverDoublets is exported by scDblFinder", {
  testthat::skip_if_not_installed("scDblFinder")
  testthat::expect_true("recoverDoublets" %in% getNamespaceExports("scDblFinder"))
})

testthat::test_that("runner forwards separate raw ADT background path", {
  runner_text <- paste(readLines(file.path(repo_root, "qc_batch_runner.R"), warn = FALSE),
                       collapse = "\n")
  testthat::expect_match(runner_text, "\"raw_adt_path\"", fixed = TRUE)
  testthat::expect_match(runner_text, "raw_adt_path     = row[[\"raw_adt_path\"]]", fixed = TRUE)
})

testthat::test_that("integration level defaults to Sample and accepts Batch", {
  runner_expressions <- parse(file.path(repo_root, "qc_batch_runner.R"))
  parser_assignment <- Filter(function(expression) {
    is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      identical(expression[[2L]], as.name("parse_cli_args"))
  }, runner_expressions)[[1L]]
  parser_environment <- new.env(parent = globalenv())
  eval(parser_assignment, envir = parser_environment)

  default_options <- parser_environment$parse_cli_args(c("--run_integration", "TRUE"))
  batch_options <- parser_environment$parse_cli_args(c(
    "--run_integration", "TRUE", "--integration_level", "Batch"
  ))
  testthat::expect_identical(default_options$integration_level, "Sample")
  testthat::expect_identical(batch_options$integration_level, "Batch")
})

testthat::test_that("shell launchers parse without submitting jobs", {
  testthat::expect_identical(system2("bash", c("-n", file.path(repo_root, "submit_QC_batch.sh"))), 0L)
  testthat::expect_identical(system2("bash", c("-n", file.path(repo_root, "run_qc_batch_local.sh"))), 0L)
})

testthat::test_that("submitter extracts sample names from transposed header", {
  command <- paste(
    paste("head -n 1", shQuote(file.path(repo_root, "samples.csv")),
          "| cut -d',' --complement -f1 |"),
    "tr ',' '\\n' | tr -d '\\r'"
  )
  sample_names <- system2("bash", c("-c", shQuote(command)), stdout = TRUE)
  testthat::expect_identical(sample_names,
                             c("Control_1", "Control_2", "Treatment_1", "Pool_1"))
})