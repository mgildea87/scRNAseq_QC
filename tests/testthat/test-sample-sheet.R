testthat::test_that("transposed sample sheet reconstructs samples and selects one", {
  sample_sheet_path <- tempfile(fileext = ".csv")
  on.exit(unlink(sample_sheet_path), add = TRUE)
  write_transposed_sample_sheet(sample_sheet_path)

  samples <- load_transposed_sample_sheet(sample_sheet_path)
  testthat::expect_identical(samples$sample_name, c("A", "B"))
  testthat::expect_identical(samples$filtered_path, c("/filtered/A", "/filtered/B"))
  testthat::expect_true(all(is.na(samples$min_nCount_RNA)))
  testthat::expect_identical(select_sample_rows(samples, "B")$sample_name, "B")
  testthat::expect_identical(nrow(select_sample_rows(samples, NULL)), 2L)
  testthat::expect_error(select_sample_rows(samples, "missing"), "not found")
})

testthat::test_that("transposed sample sheet rejects malformed structure", {
  mismatch_path <- tempfile(fileext = ".csv")
  mismatch <- data.frame(
    field = c("sample_name", "filtered_path", "raw_path"),
    A = c("Wrong", "/filtered/A", "/raw/A"),
    check.names = FALSE
  )
  write.csv(mismatch, mismatch_path, row.names = FALSE)
  on.exit(unlink(mismatch_path), add = TRUE)
  testthat::expect_error(load_transposed_sample_sheet(mismatch_path), "must match")

  duplicate_path <- tempfile(fileext = ".csv")
  duplicate <- data.frame(
    field = c("sample_name", "filtered_path", "raw_path", "raw_path"),
    A = c("A", "/filtered/A", "/raw/A", "/raw/A"),
    check.names = FALSE
  )
  write.csv(duplicate, duplicate_path, row.names = FALSE)
  on.exit(unlink(duplicate_path), add = TRUE)
  testthat::expect_error(load_transposed_sample_sheet(duplicate_path), "Duplicate field")

  missing_path <- tempfile(fileext = ".csv")
  missing <- data.frame(
    field = c("sample_name", "filtered_path"),
    A = c("A", "/filtered/A"),
    check.names = FALSE
  )
  write.csv(missing, missing_path, row.names = FALSE)
  on.exit(unlink(missing_path), add = TRUE)
  testthat::expect_error(load_transposed_sample_sheet(missing_path), "raw_path")
})