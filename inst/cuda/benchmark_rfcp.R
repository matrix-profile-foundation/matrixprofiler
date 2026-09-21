suppressPackageStartupMessages(library(matrixprofiler))

fixture_path <- "dev/mpxfixture/matrixprofiler_na_ab_fixture.rds"
if (!file.exists(fixture_path)) {
  cat("Fixture not found at", fixture_path, "\n")
  quit(status = 1)
}

fix <- readRDS(fixture_path)
pos <- fix$streams$no_na$t_plus
neg <- fix$streams$no_na$t_minus
w <- as.integer(fix$metadata$window_size)
k <- 4L

cat("========================================================================================\n")
cat("                MATRIXPROFILER RFCP BENCHMARK (dev/mpxfixture)                          \n")
cat(sprintf(" Dataset: pos = %d samples, neg = %d samples, window_size = %d, max_freq = %d\n",
            length(pos), length(neg), w, k))
cat("========================================================================================\n\n")

query_counts <- c(250L, 500L, 1000L)
results <- list()

for (qc in query_counts) {
  cat(sprintf("\n>>> Running benchmark for Q = %d queries...\n", qc))

  # 1. CPU (1 thread)
  t_cpu1 <- system.time({
    res_cpu1 <- rfcp(pos, neg, w, max_freq = k, query_begin = 1L, query_end = qc,
                     n_workers = 1L, backend = "cpu")
  })

  # 2. CPU (4 threads)
  t_cpu4 <- system.time({
    res_cpu4 <- rfcp(pos, neg, w, max_freq = k, query_begin = 1L, query_end = qc,
                     n_workers = 4L, backend = "cpu")
  })

  # 3. CPU (8 threads)
  t_cpu8 <- system.time({
    res_cpu8 <- rfcp(pos, neg, w, max_freq = k, query_begin = 1L, query_end = qc,
                     n_workers = 8L, backend = "cpu")
  })

  # 4. GPU (CUB)
  t_gpu <- system.time({
    res_gpu <- rfcp(pos, neg, w, max_freq = k, query_begin = 1L, query_end = qc,
                    backend = "cuda")
  })

  max_diff <- max(abs(res_cpu1$rms_profile - res_gpu$rms_profile), na.rm = TRUE)
  plato_match <- (res_cpu1$plato_index == res_gpu$plato_index)

  cat(sprintf("  CPU 1 thread:  %.3fs\n", t_cpu1["elapsed"]))
  cat(sprintf("  CPU 4 threads: %.3fs (speedup: %.2fx)\n", t_cpu4["elapsed"], t_cpu1["elapsed"] / t_cpu4["elapsed"]))
  cat(sprintf("  CPU 8 threads: %.3fs (speedup: %.2fx)\n", t_cpu8["elapsed"], t_cpu1["elapsed"] / t_cpu8["elapsed"]))
  cat(sprintf("  GPU (CUB):     %.3fs (speedup vs 1T: %.2fx, vs 8T: %.2fx)\n",
              t_gpu["elapsed"],
              t_cpu1["elapsed"] / t_gpu["elapsed"],
              t_cpu8["elapsed"] / t_gpu["elapsed"]))
  cat(sprintf("  Plato match:   %s (idx = %d, RMS = %.6f)\n",
              ifelse(plato_match, "YES", "NO"), res_gpu$plato_index, res_gpu$plato_rms))
  cat(sprintf("  Max abs diff:  %.3e\n", max_diff))

  results[[as.character(qc)]] <- list(
    q = qc,
    cpu1 = t_cpu1["elapsed"],
    cpu4 = t_cpu4["elapsed"],
    cpu8 = t_cpu8["elapsed"],
    gpu = t_gpu["elapsed"],
    sp_vs_1 = t_cpu1["elapsed"] / t_gpu["elapsed"],
    sp_vs_8 = t_cpu8["elapsed"] / t_gpu["elapsed"],
    max_diff = max_diff,
    plato = res_gpu$plato_index
  )
}

cat("\n\n========================================================================================\n")
cat("                               SUMMARY TABLE                                            \n")
cat("========================================================================================\n")
cat(sprintf("| %-8s | %-12s | %-12s | %-12s | %-12s | %-12s | %-12s |\n",
            "Queries", "CPU (1T)", "CPU (4T)", "CPU (8T)", "GPU (CUB)", "Sp vs 1T", "Sp vs 8T"))
cat("|----------|--------------|--------------|--------------|--------------|--------------|--------------|\n")
for (r in results) {
  cat(sprintf("| %-8d | %-10.3fs   | %-10.3fs   | %-10.3fs   | %-10.3fs   | %-10.2fx   | %-10.2fx   |\n",
              r$q, r$cpu1, r$cpu4, r$cpu8, r$gpu, r$sp_vs_1, r$sp_vs_8))
}
cat("========================================================================================\n")
