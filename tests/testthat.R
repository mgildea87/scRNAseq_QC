if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("Install the testthat package to run the local test suite.", call. = FALSE)
}

testthat::test_dir("tests/testthat", reporter = "summary")