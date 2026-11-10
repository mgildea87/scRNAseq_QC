convolution_rolling_mean <- function(values, window) {
  as.numeric(stats::filter(values, rep(1 / window, window), sides = 2))
}

make_barcode_rank_fixture <- function(counts) {
  counts <- sort(counts, decreasing = TRUE)
  n <- length(counts)
  ranks <- seq_len(n)
  log_rank <- log10(ranks)
  log_umi <- log10(pmax(counts, 1))
  window <- max(101L, as.integer(0.01 * n))
  if (window %% 2L == 0L) window <- window + 1L
  window <- min(window, n - (1L - n %% 2L))
  half_window <- (window - 1L) / 2
  left <- pmax(ranks - half_window, 1L)
  right <- pmin(ranks + half_window, n)
  slope <- rep(NA_real_, n)
  valid <- right > left
  slope[valid] <- (log_umi[right[valid]] - log_umi[left[valid]]) /
    (log_rank[right[valid]] - log_rank[left[valid]])
  data.frame(UMI_counts = counts, rank = ranks, log_rank = log_rank,
             log_umi = log_umi, moving_slope = slope)
}

testthat::test_that("linear centered mean matches convolution and keeps edge NAs", {
  set.seed(11)
  for (n in c(11L, 100L, 10001L)) {
    values <- cumsum(stats::rnorm(n)) / n
    for (window in c(1L, 11L)) {
      expected <- convolution_rolling_mean(values, window)
      actual <- centered_rolling_mean(values, window)
      testthat::expect_identical(is.na(actual), is.na(expected))
      testthat::expect_equal(actual, expected, tolerance = 1e-10)
    }
  }
  testthat::expect_equal(centered_rolling_mean(3, 1), 3)
  testthat::expect_error(centered_rolling_mean(1:10, 4), "odd integer")
  testthat::expect_error(centered_rolling_mean(c(1, NA, 2), 3), "finite")
})

testthat::test_that("shelf results match convolution for barcode curves and plateaus", {
  set.seed(4)
  counts <- c(round(10^stats::runif(10000L, 0, 4)), rep(0, 100L))
  ranked <- make_barcode_rank_fixture(counts)
  fixtures <- list(
    ranked,
    data.frame(UMI_counts = rep(10, 10000), rank = seq_len(10000),
               moving_slope = rep(-0.1, 10000)),
    data.frame(UMI_counts = rep(10, 10000), rank = seq_len(10000),
               moving_slope = -1 - 0.3 * sin(seq_len(10000) / 43))
  )

  for (data in fixtures) {
    expected <- find_ambient_shelf(data, 10, convolution_rolling_mean)
    actual <- find_ambient_shelf(data, 10)
    testthat::expect_equal(actual, expected, tolerance = 1e-10)
  }
  oscillating_result <- find_ambient_shelf(fixtures[[3L]], 10)
  testthat::expect_equal(oscillating_result$first_drop_rank, 68)
  testthat::expect_equal(oscillating_result$second_drop_rank, 338)
  testthat::expect_equal(oscillating_result$ambient_peak$rank, 203)
  testthat::expect_equal(nrow(oscillating_result$ambient_shelf), 269L)
  testthat::expect_equal(oscillating_result$average_ambient_rna_umi, 10)
})

testthat::test_that("small inputs and missing shelf boundaries use the existing fallback", {
  data <- data.frame(UMI_counts = c(9, 5, 2, 0), rank = 1:4,
                     moving_slope = c(-1, -0.5, 0, NA))
  result <- find_ambient_shelf(data, 5)
  testthat::expect_identical(result$ambient_shelf$rank, 2:3)
  testthat::expect_identical(result$ambient_peak$rank, 3L)
  testthat::expect_equal(result$average_ambient_rna_umi, 3.5)
  testthat::expect_true(is.na(result$first_drop_rank))
  testthat::expect_equal(result, find_ambient_shelf(data, 5, convolution_rolling_mean))
  testthat::expect_null(find_ambient_shelf(data, 1)$ambient_shelf)
  testthat::expect_identical(find_ambient_shelf(data, NA_real_)$ambient_shelf$rank, 1:3)
})

testthat::test_that("plot sampling is bounded, deterministic and includes landmarks", {
  testthat::expect_identical(barcode_plot_indices(3), 1:3)
  indices <- barcode_plot_indices(1000000, c(123456, 765432, NA))
  testthat::expect_identical(indices, barcode_plot_indices(1000000, c(123456, 765432, NA)))
  testthat::expect_true(all(c(1L, 1000000L, 123456L, 765432L) %in% indices))
  testthat::expect_lte(length(indices), 50002L)
  testthat::expect_true(all(diff(indices) > 0))
})

testthat::test_that("raw rank order maps tied counts to original SoupX barcodes", {
  raw_counts <- c(5, 10, 5, 2, 10)
  barcodes <- paste0("cell", seq_along(raw_counts))
  old_order <- order(raw_counts, decreasing = TRUE)
  new_order <- order(raw_counts, decreasing = TRUE, method = "radix")
  shelf_ranks <- c(2L, 3L, 4L)
  testthat::expect_identical(barcodes[new_order[shelf_ranks]],
                             barcodes[old_order[shelf_ranks]])

  template <- paste(readLines(file.path(repo_root, "sample_QC.Rmd"), warn = FALSE),
                    collapse = "\n")
  testthat::expect_match(template, "colnames(raw_matrix)[raw_umi_order[ambient_shelf$rank]]",
                         fixed = TRUE)
  testthat::expect_match(template, "find_ambient_shelf(bc_rank_data, cellranger_count_thr)",
                         fixed = TRUE)
  testthat::expect_match(template, "ggplot(bc_plot_data", fixed = TRUE)
})
