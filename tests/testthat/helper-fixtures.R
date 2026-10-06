make_synthetic_10x_fixture <- function() {
  root <- tempfile("synthetic-10x-")
  filtered_dir <- file.path(root, "filtered")
  raw_dir <- file.path(root, "raw")
  hto_only_dir <- file.path(root, "hto-only")
  adt_only_dir <- file.path(root, "adt-only")
  dir.create(filtered_dir, recursive = TRUE)
  dir.create(raw_dir, recursive = TRUE)
  dir.create(hto_only_dir, recursive = TRUE)
  dir.create(adt_only_dir, recursive = TRUE)

  feature_ids <- c("GENE_1", "GENE_2", "HTO_ID_A", "HTO_ID_B", "ADT_ID_CD3", "ADT_ID_CD4")
  feature_names <- c("GeneA", "GeneB", "HTO_A", "HTO_B", "CD3", "CD4")
  feature_types <- c(rep("Gene Expression", 2L), rep("Antibody Capture", 4L))
  features <- data.frame(id = feature_ids, name = feature_names,
                         feature_type = feature_types, stringsAsFactors = FALSE)
  barcodes <- c("cellA-1", "cellB-1", "emptyA-1", "emptyB-1")
  counts <- matrix(
    c(10, 8, 0, 0,
      5, 3, 0, 0,
      4, 3, 1, 0,
      0, 0, 1, 0,
      3, 1, 1, 2,
      2, 4, 0, 1),
    nrow = length(feature_names),
    byrow = TRUE,
    dimnames = list(feature_names, barcodes)
  )

  write_matrix <- function(directory, matrix_counts, matrix_barcodes,
                           matrix_features = features) {
    plain_paths <- file.path(directory, c("matrix.mtx", "barcodes.tsv", "features.tsv"))
    Matrix::writeMM(Matrix::Matrix(matrix_counts, sparse = TRUE), plain_paths[[1L]])
    write.table(matrix_barcodes, plain_paths[[2L]],
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(matrix_features, plain_paths[[3L]],
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

    for (plain_path in plain_paths) {
      compressed_path <- paste0(plain_path, ".gz")
      compressed <- gzfile(compressed_path, open = "wt")
      writeLines(readLines(plain_path, warn = FALSE), compressed)
      close(compressed)
      unlink(plain_path)
    }
  }

  write_matrix(filtered_dir, counts[, 1:2, drop = FALSE], barcodes[1:2])
  write_matrix(raw_dir, counts, barcodes)
  hto_features_table <- features[3:4, , drop = FALSE]
  hto_features_table$feature_type <- "Multiplexing Capture"
  write_matrix(hto_only_dir, counts[3:4, 1:2, drop = FALSE], barcodes[1:2],
               hto_features_table)
  write_matrix(adt_only_dir, counts[5:6, 1:2, drop = FALSE], barcodes[1:2],
               features[5:6, , drop = FALSE])

  feature_reference_path <- file.path(root, "feature_reference.csv")
  write.csv(features[3:6, ], feature_reference_path, row.names = FALSE)

  list(
    root = root,
    filtered_dir = filtered_dir,
    raw_dir = raw_dir,
    hto_only_dir = hto_only_dir,
    adt_only_dir = adt_only_dir,
    feature_reference_path = feature_reference_path,
    feature_names = feature_names[3:6],
    hto_features = "HTO_A;HTO_B"
  )
}

write_transposed_sample_sheet <- function(path, sample_names = c("A", "B")) {
  sheet <- data.frame(
    field = c("sample_name", "filtered_path", "raw_path", "batch",
              "min_nCount_RNA", "use_hashtag"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  for (sample_name in sample_names) {
    sheet[[sample_name]] <- c(sample_name, paste0("/filtered/", sample_name),
                              paste0("/raw/", sample_name), "batch1", NA, "FALSE")
  }
  write.csv(sheet, path, row.names = FALSE, na = "NA", quote = TRUE)
}