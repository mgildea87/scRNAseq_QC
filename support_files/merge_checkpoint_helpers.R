atomic_save_rds <- function(object, path) {
  temporary_path <- tempfile(pattern = paste0(".", basename(path), "-"),
                             tmpdir = dirname(path))
  on.exit(unlink(temporary_path), add = TRUE)
  saveRDS(object, temporary_path)
  if (!file.rename(temporary_path, path)) {
    stop("Could not replace checkpoint: ", path, call. = FALSE)
  }
}

checkpoint_file_hash <- function(path) {
  hash <- unname(tools::md5sum(path))
  if (length(hash) != 1L || is.na(hash)) {
    stop("Could not checksum checkpoint: ", path, call. = FALSE)
  }
  hash
}

checkpoint_fingerprint <- function(inputs) {
  temporary_path <- tempfile("lisi-inputs-")
  on.exit(unlink(temporary_path), add = TRUE)
  saveRDS(inputs, temporary_path, compress = FALSE)
  checkpoint_file_hash(temporary_path)
}

rna_layers_need_join <- function(assay) {
  any(vapply(c("counts", "data", "scale.data"), function(layer) {
    matches <- suppressWarnings(SeuratObject::Layers(assay, search = layer))
    length(matches) > 1L || (length(matches) == 1L && matches[[1L]] != layer)
  }, logical(1)))
}

get_lisi_checkpoint <- function(path, fingerprint, reuse, compute_scores) {
  provenance_path <- paste0(path, ".provenance.rds")
  if (reuse && file.exists(path) && file.exists(provenance_path)) {
    provenance <- readRDS(provenance_path)
    if (!is.list(provenance) ||
        !identical(provenance$version, 1L) ||
        !is.character(provenance$input_hash) || length(provenance$input_hash) != 1L ||
        is.na(provenance$input_hash) ||
        !is.character(provenance$result_hash) || length(provenance$result_hash) != 1L ||
        is.na(provenance$result_hash)) {
      stop("Invalid LISI checkpoint provenance: ", provenance_path, call. = FALSE)
    }
    if (identical(provenance$input_hash, fingerprint) &&
        identical(provenance$result_hash, checkpoint_file_hash(path))) {
      message("Reusing verified LISI checkpoint: ", path)
      return(readRDS(path))
    }
    message("LISI checkpoint inputs or result do not match; recalculating.")
  } else if (reuse) {
    message("LISI checkpoint or provenance is missing; recalculating.")
  }

  result <- compute_scores()
  atomic_save_rds(result, path)
  atomic_save_rds(list(
    version = 1L,
    input_hash = fingerprint,
    result_hash = checkpoint_file_hash(path)
  ), provenance_path)
  result
}
