suppressPackageStartupMessages({
  library(matrixprofiler)
})

cat("========================================================\n")
cat("          MATRIXPROFILER RFCP GPU VERIFICATION          \n")
cat("========================================================\n\n")

# ---------------------------------------------------------------------------
# Test 1: Periodic series with NA gap (basic validation & NA robustness)
# ---------------------------------------------------------------------------
cat("Test 1: Periodic series with NA gap...\n")
set.seed(42)
pos_sine <- c(sin(seq(0, 4 * pi, length.out = 200)), NA, sin(seq(0, 4 * pi, length.out = 200)))
neg_sine <- cos(seq(0, 10 * pi, length.out = 450))
m1 <- 20L
k1 <- 4L

t_cpu1 <- system.time(cpu1 <- rfcp(pos_sine, neg_sine, m1, max_freq = k1, backend = "cpu"))
t_gpu1 <- system.time(gpu1 <- rfcp(pos_sine, neg_sine, m1, max_freq = k1, backend = "cuda"))

cat(sprintf("  Time: CPU = %.3fs, GPU = %.3fs\n", t_cpu1["elapsed"], t_gpu1["elapsed"]))
cat(sprintf("  Plato index: CPU = %d, GPU = %d\n", cpu1$plato_index, gpu1$plato_index))
cat(sprintf("  Plato RMS:   CPU = %.8f, GPU = %.8f\n", cpu1$plato_rms, gpu1$plato_rms))

diff_rms1 <- max(abs(cpu1$rms_profile - gpu1$rms_profile), na.rm = TRUE)
cat(sprintf("  Max abs diff in rms_profile: %.3e\n", diff_rms1))

stopifnot(cpu1$plato_index == gpu1$plato_index)
stopifnot(abs(cpu1$plato_rms - gpu1$plato_rms) < 1e-6)
stopifnot(diff_rms1 < 1e-6)
stopifnot(identical(which(is.na(cpu1$rms_profile)), which(is.na(gpu1$rms_profile))))
cat("  -> Test 1 PASSED!\n\n")

# ---------------------------------------------------------------------------
# Test 2: Realistic noisy time series with return_profiles = TRUE
# ---------------------------------------------------------------------------
cat("Test 2: Realistic random walk series with return_profiles = TRUE...\n")
set.seed(123)
pos_rw <- cumsum(rnorm(300))
pos_rw[150] <- NA
neg_rw <- cumsum(rnorm(350))
m2 <- 25L
k2 <- 4L

t_cpu2 <- system.time(
  cpu2 <- rfcp(pos_rw, neg_rw, m2, max_freq = k2, return_profiles = TRUE, backend = "cpu")
)
t_gpu2 <- system.time(
  gpu2 <- rfcp(pos_rw, neg_rw, m2, max_freq = k2, return_profiles = TRUE, backend = "cuda")
)

cat(sprintf("  Time: CPU = %.3fs, GPU = %.3fs\n", t_cpu2["elapsed"], t_gpu2["elapsed"]))
cat(sprintf("  Plato index: CPU = %d, GPU = %d\n", cpu2$plato_index, gpu2$plato_index))
cat(sprintf("  Plato RMS:   CPU = %.8f, GPU = %.8f\n", cpu2$plato_rms, gpu2$plato_rms))

diff_rms2 <- max(abs(cpu2$rms_profile - gpu2$rms_profile), na.rm = TRUE)
diff_rfmp_aa <- max(abs(cpu2$rfmp_aa - gpu2$rfmp_aa), na.rm = TRUE)
diff_rfmp_ab <- max(abs(cpu2$rfmp_ab - gpu2$rfmp_ab), na.rm = TRUE)
diff_rfcp <- max(abs(cpu2$rfcp - gpu2$rfcp), na.rm = TRUE)

mismatched_aa <- sum(cpu2$rfmpi_aa != gpu2$rfmpi_aa, na.rm = TRUE)
mismatched_ab <- sum(cpu2$rfmpi_ab != gpu2$rfmpi_ab, na.rm = TRUE)

cat(sprintf("  Max abs diff rms_profile: %.3e\n", diff_rms2))
cat(sprintf("  Max abs diff rfmp_aa:     %.3e\n", diff_rfmp_aa))
cat(sprintf("  Max abs diff rfmp_ab:     %.3e\n", diff_rfmp_ab))
cat(sprintf("  Max abs diff rfcp:        %.3e\n", diff_rfcp))
cat(sprintf("  Mismatched rfmpi_aa cells: %d / %d\n", mismatched_aa, length(cpu2$rfmpi_aa)))
cat(sprintf("  Mismatched rfmpi_ab cells: %d / %d\n", mismatched_ab, length(cpu2$rfmpi_ab)))

stopifnot(cpu2$plato_index == gpu2$plato_index)
stopifnot(diff_rms2 < 1e-10)
stopifnot(diff_rfmp_aa < 1e-10)
stopifnot(diff_rfmp_ab < 1e-10)
stopifnot(diff_rfcp < 1e-10)
stopifnot(mismatched_aa == 0)
stopifnot(mismatched_ab == 0)
cat("  -> Test 2 PASSED!\n\n")

# ---------------------------------------------------------------------------
# Test 3: Query subsetting (query_begin and query_end)
# ---------------------------------------------------------------------------
cat("Test 3: Query subsetting (query_begin = 50, query_end = 120)...\n")
q_begin <- 50L
q_end <- 120L
cpu3 <- rfcp(pos_rw, neg_rw, m2, max_freq = k2, query_begin = q_begin, query_end = q_end, backend = "cpu")
gpu3 <- rfcp(pos_rw, neg_rw, m2, max_freq = k2, query_begin = q_begin, query_end = q_end, backend = "cuda")

diff_subset <- max(abs(cpu3$rms_profile - gpu3$rms_profile), na.rm = TRUE)
cat(sprintf("  Max abs diff in subset rms_profile: %.3e\n", diff_subset))
stopifnot(diff_subset < 1e-10)
stopifnot(length(gpu3$rms_profile) == (q_end - q_begin + 1L))
cat("  -> Test 3 PASSED!\n\n")

cat("========================================================\n")
cat("           ALL TESTS COMPLETED SUCCESSFULLY!            \n")
cat("========================================================\n")
