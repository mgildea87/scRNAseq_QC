root_candidates <- normalizePath(c(".", "..", "../.."), mustWork = FALSE)
repo_root <- root_candidates[file.exists(file.path(root_candidates, "support_files", "qc_helpers.R"))][[1L]]
sys.source(file.path(repo_root, "support_files", "qc_helpers.R"), envir = environment())