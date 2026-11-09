testthat::test_that("SoupX contamination accepts automatic and manual rho", {
  testthat::expect_true(is.na(resolve_soupx_contamination()))
  testthat::expect_true(is.na(resolve_soupx_contamination(NA_character_)))
  testthat::expect_true(is.na(resolve_soupx_contamination("NA")))
  testthat::expect_equal(resolve_soupx_contamination(0.08), 0.08)
  testthat::expect_equal(resolve_soupx_contamination("0.08"), 0.08)
  testthat::expect_error(resolve_soupx_contamination(1), "\\[0, 1\\)")
  testthat::expect_error(resolve_soupx_contamination(-0.1), "\\[0, 1\\)")
  testthat::expect_error(resolve_soupx_contamination("invalid"), "\\[0, 1\\)")
})

testthat::test_that("SoupX rho values are resolved independently per sample", {
  rho <- normalize_sample_soupx_contamination(
    c("NA", NA_character_, "0.12", "", "0.35"),
    c("automatic", "missing", "low_rho", "blank_defaults", "high_rho")
  )

  testthat::expect_true(is.na(rho[[1L]]))
  testthat::expect_true(is.na(rho[[2L]]))
  testthat::expect_equal(rho[[3L]], 0.12)
  testthat::expect_true(is.na(rho[[4L]]))
  testthat::expect_equal(rho[[5L]], 0.35)
  testthat::expect_error(
    normalize_sample_soupx_contamination(c("NA", "1"), c("sample_a", "sample_b")),
    "sample 'sample_b'"
  )
})

testthat::test_that("SoupX matrices share features and retain raw background droplets", {
  raw_counts <- matrix(
    1:12, nrow = 3,
    dimnames = list(c("gene_a", "gene_b", "gene_c"), c("cell_1", "cell_2", "empty_1", "empty_2"))
  )
  filtered_counts <- matrix(
    1:4, nrow = 2,
    dimnames = list(c("gene_b", "gene_a"), c("cell_1", "cell_2"))
  )

  aligned <- align_soupx_matrices(raw_counts, filtered_counts)
  testthat::expect_identical(rownames(aligned$tod), c("gene_b", "gene_a"))
  testthat::expect_identical(colnames(aligned$tod), colnames(raw_counts))
  testthat::expect_identical(rownames(aligned$toc), c("gene_b", "gene_a"))
  testthat::expect_identical(colnames(aligned$toc), colnames(filtered_counts))
  filtered_with_unknown_cell <- cbind(
    filtered_counts,
    matrix(c(5, 6), ncol = 1, dimnames = list(rownames(filtered_counts), "cell_3"))
  )
  testthat::expect_error(
    align_soupx_matrices(raw_counts, filtered_with_unknown_cell),
    "all filtered barcodes"
  )
})
