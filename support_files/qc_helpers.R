read_feature_reference_names <- function(feature_reference_path) {
  if (length(feature_reference_path) != 1L || is.na(feature_reference_path) ||
      !nzchar(trimws(feature_reference_path)) || !file.exists(feature_reference_path)) {
    stop("feature_reference_path must point to an existing file.", call. = FALSE)
  }

  feature_reference <- read.csv(feature_reference_path, stringsAsFactors = FALSE)
  if (!"name" %in% colnames(feature_reference)) {
    stop("Feature reference is missing the required 'name' column: ",
         feature_reference_path, call. = FALSE)
  }

  feature_names <- unique(trimws(as.character(feature_reference$name)))
  feature_names <- feature_names[!is.na(feature_names) & nzchar(feature_names)]
  if (length(feature_names) == 0L) {
    stop("Feature reference contains no non-empty feature names: ",
         feature_reference_path, call. = FALSE)
  }

  feature_names
}

select_10x_assay_matrix <- function(x, feature_types, label, input_path) {
  if (inherits(x, "Matrix") || is.matrix(x)) return(x)
  if (!is.list(x) || length(x) == 0L) {
    stop("Unsupported ", label, " input format at: ", input_path, call. = FALSE)
  }
  if (length(x) == 1L) return(x[[1L]])
  if (is.null(names(x)) || !any(nzchar(names(x)))) {
    stop("Multiple unnamed matrices detected in ", label, " input: ", input_path,
         call. = FALSE)
  }

  for (feature_type in feature_types) {
    match_idx <- which(tolower(names(x)) == tolower(feature_type))
    if (length(match_idx) > 0L) return(x[[match_idx[[1L]]]])
  }

  stop(
    "Could not find a requested feature type (", paste(feature_types, collapse = ", "),
    ") in ", label, " input: ", input_path,
    call. = FALSE
  )
}

validate_raw_adt_path <- function(use_adt, raw_adt_path) {
  if (!isTRUE(use_adt)) return(invisible(NULL))

  if (length(raw_adt_path) != 1L) {
    stop("raw_adt_path must be a single path when provided.", call. = FALSE)
  }
  if (is.na(raw_adt_path) || !nzchar(trimws(raw_adt_path))) return(invisible(NULL))
  if (!file.exists(raw_adt_path)) {
    stop("raw_adt_path must point to an existing raw ADT matrix when provided.",
         call. = FALSE)
  }

  invisible(raw_adt_path)
}

split_antibody_capture <- function(ab_matrix, feature_names, hto_features) {
  if (is.null(rownames(ab_matrix))) {
    stop("Antibody Capture matrix must have feature row names.", call. = FALSE)
  }

  missing_matrix_names <- setdiff(feature_names, rownames(ab_matrix))
  if (length(missing_matrix_names) > 0L) {
    stop("Feature reference names missing from Antibody Capture matrix: ",
         paste(missing_matrix_names, collapse = ", "), call. = FALSE)
  }

  explicit_hto_names <- trimws(strsplit(as.character(hto_features), "[;,]")[[1L]])
  explicit_hto_names <- unique(explicit_hto_names[nzchar(explicit_hto_names)])
  missing_hto_names <- setdiff(explicit_hto_names, feature_names)
  if (length(missing_hto_names) > 0L) {
    stop("HTO feature names not found in feature reference: ",
         paste(missing_hto_names, collapse = ", "), call. = FALSE)
  }

  hto_names <- intersect(explicit_hto_names, feature_names)
  adt_names <- setdiff(feature_names, hto_names)

  list(
    hto = if (length(hto_names) > 0L) ab_matrix[hto_names, , drop = FALSE] else NULL,
    adt = if (length(adt_names) > 0L) ab_matrix[adt_names, , drop = FALSE] else NULL
  )
}

build_mad_table <- function(values, med, mad) {
  mad_level <- 5:-5
  empirical_cdf <- ecdf(values)
  filter_threshold <- med + mad_level * mad
  percentiles <- empirical_cdf(filter_threshold)
  data.frame(
    mad = mad_level,
    filter_threshold = filter_threshold,
    percent_below = percentiles,
    cells_below = percentiles * length(values),
    percent_above = 1 - percentiles,
    cells_above = (1 - percentiles) * length(values)
  )
}

resolve_adt_count_thresholds <- function(values, min_override = NA_real_, max_override = NA_real_) {
  values <- values[!is.na(values)]
  if (length(values) == 0L) {
    stop("Cannot resolve ADT count thresholds without non-missing values.", call. = FALSE)
  }

  median_value <- stats::median(values)
  mad_value <- stats::median(abs(values - median_value))
  list(
    median = median_value,
    mad = mad_value,
    min = if (is.na(min_override)) median_value - 4 * mad_value else min_override,
    max = if (is.na(max_override)) Inf else max_override
  )
}

adt_count_keep_mask <- function(values, min_threshold, max_threshold, enabled = TRUE) {
  if (!isTRUE(enabled)) return(rep(TRUE, length(values)))
  !is.na(values) & values > min_threshold & values <= max_threshold
}

hto_singlet_keep_mask <- function(hto_class) {
  !is.na(hto_class) & as.character(hto_class) == "Singlet"
}

resolve_hto_doublet_dimreds <- function(reduced_dim_names) {
  reduced_dim_names <- as.character(reduced_dim_names)
  pca_names <- reduced_dim_names[grepl("pca", reduced_dim_names, ignore.case = TRUE)]
  umap_names <- reduced_dim_names[grepl("umap", reduced_dim_names, ignore.case = TRUE)]
  if (length(pca_names) == 0L || length(umap_names) == 0L) {
    stop("SingleCellExperiment conversion must preserve PCA and UMAP reduced dimensions.",
         call. = FALSE)
  }

  list(pca = pca_names[[1L]], umap = umap_names[[1L]])
}

hto_singlet_expression_thresholds <- function(expression_matrix, hto_class, hash_id) {
  feature_names <- rownames(expression_matrix)
  cell_names <- colnames(expression_matrix)
  if (is.null(feature_names) || is.null(cell_names)) {
    stop("HTO expression matrix must have feature and cell names.", call. = FALSE)
  }

  if (!is.null(names(hto_class)) && !is.null(names(hash_id))) {
    hto_class <- hto_class[cell_names]
    hash_id <- hash_id[cell_names]
  }
  if (length(hto_class) != length(cell_names) || length(hash_id) != length(cell_names)) {
    stop("HTO classification metadata must match the expression matrix cells.", call. = FALSE)
  }

  hto_class <- as.character(hto_class)
  hash_id <- as.character(hash_id)
  thresholds <- stats::setNames(rep(NA_real_, length(feature_names)), feature_names)
  for (feature_name in feature_names) {
    singlet_cells <- !is.na(hto_class) & hto_class == "Singlet" &
      !is.na(hash_id) & hash_id == feature_name
    expression_values <- as.numeric(expression_matrix[feature_name, singlet_cells, drop = TRUE])
    expression_values <- expression_values[is.finite(expression_values)]
    if (length(expression_values) > 0L) {
      thresholds[[feature_name]] <- min(expression_values)
    }
  }

  thresholds
}

load_transposed_sample_sheet <- function(path,
                                         required_fields = c("sample_name", "filtered_path", "raw_path")) {
  sheet <- read.csv(
    path,
    stringsAsFactors = FALSE,
    strip.white = TRUE,
    check.names = FALSE,
    colClasses = "character"
  )
  if (ncol(sheet) < 2L || !identical(names(sheet)[[1L]], "field")) {
    stop("Transposed sample sheet must start with a 'field' column and at least one sample column.",
         call. = FALSE)
  }

  fields <- trimws(sheet[[1L]])
  if (anyNA(fields) || any(!nzchar(fields))) {
    stop("Sample sheet contains an empty field name.", call. = FALSE)
  }
  if (anyDuplicated(fields)) {
    stop("Duplicate field names found in transposed sample sheet.", call. = FALSE)
  }

  sample_names <- names(sheet)[-1L]
  if (anyNA(sample_names) || any(!nzchar(sample_names)) || anyDuplicated(sample_names)) {
    stop("Sample columns must have non-empty, unique names.", call. = FALSE)
  }

  missing_fields <- setdiff(required_fields, fields)
  if (length(missing_fields) > 0L) {
    stop("Sample sheet is missing required field(s): ",
         paste(missing_fields, collapse = ", "), call. = FALSE)
  }

  sample_sheet <- as.data.frame(t(as.matrix(sheet[-1L])), stringsAsFactors = FALSE,
                                optional = TRUE)
  names(sample_sheet) <- fields
  rownames(sample_sheet) <- NULL

  if (!identical(as.character(sample_sheet$sample_name), sample_names)) {
    stop("Each sample_name field value must match its sample column header.", call. = FALSE)
  }
  for (field in required_fields) {
    values <- sample_sheet[[field]]
    if (anyNA(values) || any(!nzchar(trimws(values)))) {
      stop("Required field '", field, "' contains a blank value.", call. = FALSE)
    }
  }
  if (anyDuplicated(sample_sheet$sample_name)) {
    stop("Duplicate sample_name entries found in sample sheet.", call. = FALSE)
  }

  sample_sheet
}

select_sample_rows <- function(sample_sheet, sample_name) {
  if (is.null(sample_name)) return(sample_sheet)

  selected <- sample_sheet[sample_sheet$sample_name == sample_name, , drop = FALSE]
  if (nrow(selected) == 0L) {
    stop("Sample '", sample_name, "' not found in sample sheet.", call. = FALSE)
  }
  selected
}