testthat::test_that("synthetic CellRanger MEX inputs pass through Read10X", {
  testthat::skip_if_not_installed("Seurat")
  fixture <- make_synthetic_10x_fixture()
  on.exit(unlink(fixture$root, recursive = TRUE), add = TRUE)

  filtered <- Seurat::Read10X(fixture$filtered_dir)
  raw <- Seurat::Read10X(fixture$raw_dir)
  testthat::expect_true(is.list(filtered))
  testthat::expect_true(all(c("Gene Expression", "Antibody Capture") %in% names(filtered)))
  testthat::expect_identical(rownames(filtered[["Gene Expression"]]), c("GeneA", "GeneB"))

  reference_names <- read_feature_reference_names(fixture$feature_reference_path)
  split_filtered <- split_antibody_capture(
    filtered[["Antibody Capture"]], reference_names, fixture$hto_features
  )
  split_raw <- split_antibody_capture(
    raw[["Antibody Capture"]], reference_names, fixture$hto_features
  )
  testthat::expect_identical(rownames(split_filtered$hto), c("HTO_A", "HTO_B"))
  testthat::expect_identical(rownames(split_filtered$adt), c("CD3", "CD4"))
  testthat::expect_identical(rownames(split_raw$adt), c("CD3", "CD4"))
  testthat::expect_true(all(c("emptyA-1", "emptyB-1") %in%
                              setdiff(colnames(raw[["Gene Expression"]]),
                                      colnames(filtered[["Gene Expression"]]))))

  hto_only <- select_10x_assay_matrix(
    Seurat::Read10X(fixture$hto_only_dir),
    c("Multiplexing Capture", "Antibody Capture"),
    "HTO", fixture$hto_only_dir
  )
  adt_only <- select_10x_assay_matrix(
    Seurat::Read10X(fixture$adt_only_dir),
    "Antibody Capture",
    "ADT", fixture$adt_only_dir
  )
  testthat::expect_identical(rownames(hto_only), c("HTO_A", "HTO_B"))
  testthat::expect_identical(rownames(adt_only), c("CD3", "CD4"))
})