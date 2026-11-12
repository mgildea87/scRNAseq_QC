testthat::test_that("raw matrix HTO names match Seurat-normalized assignments", {
  raw_names <- paste0("Totalseq_A0", 301:310)
  raw <- matrix(1, nrow = length(raw_names), ncol = 2,
                dimnames = list(raw_names, c("cell1", "cell2")))
  assay_names <- rownames(dash_feature_names(raw))
  calls <- c("Totalseq-A0301", "Totalseq-A0305", "Totalseq-A0306", NA)
  keep <- paste(raw_names[1:5], collapse = ";")

  testthat::expect_identical(
    hto_id_keep_mask(calls, keep, assay_names),
    c(TRUE, TRUE, FALSE, FALSE)
  )
  testthat::expect_identical(
    hto_id_keep_mask(calls, "Totalseq-A0301, Totalseq_A0305", assay_names),
    c(TRUE, TRUE, FALSE, FALSE)
  )
  testthat::expect_identical(hto_id_keep_mask(calls, "", assay_names), rep(TRUE, 4))
  testthat::expect_error(
    hto_id_keep_mask(character(), "Totalseq_A0999", assay_names),
    "Totalseq_A0999"
  )
})

testthat::test_that("HTO keep-list is validated before doublet detection", {
  template <- paste(readLines(file.path(repo_root, "sample_QC.Rmd"), warn = FALSE),
                    collapse = "\n")
  validation <- regexpr("hto_id_keep_mask(character(), params$hto_ids_to_keep", template,
                        fixed = TRUE)[[1L]]
  doublets <- regexpr("```{r hto-intra-sample-doublets", template, fixed = TRUE)[[1L]]
  filtering <- regexpr("```{r hto-filter-selected-ids", template, fixed = TRUE)[[1L]]
  testthat::expect_gt(validation, 0L)
  testthat::expect_lt(validation, doublets)
  testthat::expect_lt(doublets, filtering)
})
