testthat::test_that("ADT thresholds use MAD defaults and honor explicit bounds", {
  values <- c(0, 10, 20, 30, 40)
  defaults <- resolve_adt_count_thresholds(values)

  testthat::expect_equal(defaults$median, 20)
  testthat::expect_equal(defaults$mad, 10)
  testthat::expect_equal(defaults$min, -20)
  testthat::expect_true(is.infinite(defaults$max))
  testthat::expect_equal(nrow(build_mad_table(values, defaults$median, defaults$mad)), 11L)

  explicit <- resolve_adt_count_thresholds(values, min_override = 10, max_override = 30)
  testthat::expect_identical(
    adt_count_keep_mask(values, explicit$min, explicit$max),
    c(FALSE, FALSE, TRUE, TRUE, FALSE)
  )
})

testthat::test_that("ADT filtering is a no-op when not enabled", {
  values <- c(NA, 0, 10, 20)
  testthat::expect_identical(
    adt_count_keep_mask(values, min_threshold = 10, max_threshold = 20, enabled = FALSE),
    rep(TRUE, length(values))
  )
})

testthat::test_that("HTO keep mask retains Singlets only", {
  testthat::expect_identical(
    hto_singlet_keep_mask(c("Singlet", "Doublet", "Negative", NA)),
    c(TRUE, FALSE, FALSE, FALSE)
  )
})

testthat::test_that("HTO doublet check resolves converted PCA and UMAP names", {
  testthat::expect_identical(
    resolve_hto_doublet_dimreds(c("PCA", "UMAP")),
    list(pca = "PCA", umap = "UMAP")
  )
  testthat::expect_identical(
    resolve_hto_doublet_dimreds(c("pca_hto_check", "umap_hto_check")),
    list(pca = "pca_hto_check", umap = "umap_hto_check")
  )
  testthat::expect_error(resolve_hto_doublet_dimreds("UMAP"), "must preserve PCA and UMAP")
})

testthat::test_that("raw ADT path is an optional override when ADT is enabled", {
  validate_raw_adt_path(FALSE, "")
  validate_raw_adt_path(TRUE, "")
  validate_raw_adt_path(TRUE, NA_character_)
  testthat::expect_error(validate_raw_adt_path(TRUE, "/missing/raw_adt"), "raw_adt_path")

  raw_adt_path <- tempfile()
  file.create(raw_adt_path)
  on.exit(unlink(raw_adt_path), add = TRUE)
  validate_raw_adt_path(TRUE, raw_adt_path)
})
testthat::test_that("dash_feature_names matches Seurat renaming so raw and filtered ADT index alike", {
  testthat::skip_if_not_installed("Seurat")
  raw <- Matrix::Matrix(
    matrix(1:12, nrow = 3, dimnames = list(c("CD3_TotalSeqB", "CD19_TotalSeqB", "Isotype_Ctrl"), paste0("bc", 1:4))),
    sparse = TRUE
  )
  renamed <- dash_feature_names(raw)
  assay <- suppressWarnings(SeuratObject::CreateAssayObject(counts = renamed[, 1:2]))

  testthat::expect_identical(rownames(assay), rownames(renamed))
  testthat::expect_equal(unname(as.matrix(renamed[rownames(assay), 3:4])), unname(as.matrix(raw[, 3:4])))
  testthat::expect_null(dash_feature_names(NULL))
  testthat::expect_error(raw[rownames(assay), , drop = FALSE], "subscript out of bounds")
})

testthat::test_that("dash_feature_names rejects colliding names", {
  m <- matrix(1:4, nrow = 2, dimnames = list(c("A_B", "A-B"), c("x", "y")))
  testthat::expect_error(dash_feature_names(m), "collide")
})
