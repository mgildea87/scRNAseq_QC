testthat::test_that("HTO keep-list defaults to all IDs and subsets when supplied", {
  hash_ids <- c("HTO_A", "HTO_B", "HTO_A", "HTO_C", NA)
  available_ids <- c("HTO_A", "HTO_B", "HTO_C")

  testthat::expect_identical(
    hto_id_keep_mask(hash_ids, "", available_ids),
    c(TRUE, TRUE, TRUE, TRUE, TRUE)
  )
  testthat::expect_identical(
    hto_id_keep_mask(hash_ids, "HTO_A; HTO_C", available_ids),
    c(TRUE, FALSE, TRUE, TRUE, FALSE)
  )
  testthat::expect_identical(
    hto_id_keep_mask(hash_ids, "HTO_B,HTO_C", available_ids),
    c(FALSE, TRUE, FALSE, TRUE, FALSE)
  )
})

testthat::test_that("HTO keep-list rejects malformed and unknown IDs", {
  testthat::expect_error(
    hto_id_keep_mask("HTO_A", " , ; ", "HTO_A"),
    "must contain at least one HTO ID"
  )
  testthat::expect_error(
    hto_id_keep_mask("HTO_A", "HTO_UNKNOWN", "HTO_A"),
    "not found in the HTO assay"
  )
  testthat::expect_error(
    hto_id_keep_mask("HTO_A", c("HTO_A", "HTO_B"), "HTO_A"),
    "single comma- or semicolon-separated string"
  )
})

testthat::test_that("runner forwards hto_ids_to_keep and sample sheet includes it", {
  runner_text <- paste(
    readLines(file.path(repo_root, "qc_batch_runner.R"), warn = FALSE),
    collapse = "\n"
  )
  template_text <- paste(
    readLines(file.path(repo_root, "sample_QC.Rmd"), warn = FALSE),
    collapse = "\n"
  )
  sample_sheet <- read.csv(
    file.path(repo_root, "samples.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  testthat::expect_true("hto_ids_to_keep" %in% sample_sheet$field)
  default_keep_list <- sample_sheet$Pool_1[sample_sheet$field == "hto_ids_to_keep"]
  testthat::expect_true(is.na(default_keep_list) || !nzchar(default_keep_list))
  testthat::expect_match(runner_text, '"hto_ids_to_keep"', fixed = TRUE)
  testthat::expect_match(runner_text, 'hto_ids_to_keep  = row[["hto_ids_to_keep"]]', fixed = TRUE)
  testthat::expect_match(template_text, "hto_ids_to_keep = params$hto_ids_to_keep", fixed = TRUE)
  singlet_filter_position <- regexpr("hto-filter-singlets", template_text, fixed = TRUE)[[1L]]
  selected_filter_position <- regexpr("hto-filter-selected-ids", template_text, fixed = TRUE)[[1L]]
  testthat::expect_gt(singlet_filter_position, 0L)
  testthat::expect_gt(selected_filter_position, singlet_filter_position)
})
