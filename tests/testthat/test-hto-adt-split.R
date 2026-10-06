testthat::test_that("antibody features split by canonical feature-reference names", {
  fixture <- make_synthetic_10x_fixture()
  on.exit(unlink(fixture$root, recursive = TRUE), add = TRUE)

  feature_names <- read_feature_reference_names(fixture$feature_reference_path)
  antibody_matrix <- matrix(
    seq_len(length(feature_names) * 2L),
    nrow = length(feature_names),
    dimnames = list(feature_names, c("cellA", "cellB"))
  )
  split <- split_antibody_capture(antibody_matrix, feature_names, fixture$hto_features)

  testthat::expect_identical(rownames(split$hto), c("HTO_A", "HTO_B"))
  testthat::expect_identical(rownames(split$adt), c("CD3", "CD4"))
  testthat::expect_equal(sum(split$hto), sum(antibody_matrix[c("HTO_A", "HTO_B"), ]))
  testthat::expect_equal(sum(split$adt), sum(antibody_matrix[c("CD3", "CD4"), ]))
})

testthat::test_that("invalid feature-reference and HTO names fail clearly", {
  test_matrix <- matrix(1:3, nrow = 3,
                        dimnames = list(c("HTO_A", "CD3", "CD4"), "cellA"))
  testthat::expect_error(
    split_antibody_capture(test_matrix, c("HTO_A", "CD3", "CD4", "CD8"), "HTO_A"),
    "missing from Antibody Capture matrix"
  )
  testthat::expect_error(
    split_antibody_capture(test_matrix, c("HTO_A", "CD3", "CD4"), "HTO_UNKNOWN"),
    "not found in feature reference"
  )

  invalid_reference <- tempfile(fileext = ".csv")
  write.csv(data.frame(id = "HTO_ID_A"), invalid_reference, row.names = FALSE)
  on.exit(unlink(invalid_reference), add = TRUE)
  testthat::expect_error(read_feature_reference_names(invalid_reference), "required 'name' column")
})

testthat::test_that("HTO thresholds use minimum CLR expression among assigned Singlets", {
  expression <- rbind(
    HTO_A = c(0.2, 0.8, 0.5, 0.1, 0.3),
    HTO_B = c(0.1, 0.2, 0.6, 0.9, 0.4),
    HTO_C = c(0.1, 0.2, 0.3, 0.4, 0.5)
  )
  colnames(expression) <- paste0("cell", seq_len(ncol(expression)))
  thresholds <- hto_singlet_expression_thresholds(
    expression,
    hto_class = c("Singlet", "Singlet", "Singlet", "Doublet", "Singlet"),
    hash_id = c("HTO_A", "HTO_A", "HTO_B", "HTO_A", "Other")
  )

  testthat::expect_equal(thresholds[["HTO_A"]], 0.2)
  testthat::expect_equal(thresholds[["HTO_B"]], 0.6)
  testthat::expect_true(is.na(thresholds[["HTO_C"]]))
})