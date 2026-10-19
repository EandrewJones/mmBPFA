# Sweep-time scaling of mmBPFA_sampler() across N, P, sparsity and threads.
# Usage: Rscript bench/scaling.R [output.csv]
# Data are drawn from the model (K = 5, binary) directly in long format, so
# sparse settings never build the N x P matrix. Time per sweep is measured
# inside the chain, so data preparation and starting values are excluded.
suppressMessages(library(mmBPFA))

simulate_long <- function(n, p, sparsity, K = 5, seed = 1) {
  set.seed(seed)
  nnz <- round((1 - sparsity) * n * p)
  cells <- sample.int(n * p, nnz)
  row <- (cells - 1) %% n + 1; col <- (cells - 1) %/% n + 1
  L <- matrix(rnorm(p * K), p, K) * matrix(rbinom(p * K, 1, 0.4), p, K)
  O <- matrix(rnorm(n * K), n, K); a <- rnorm(p, 0, 0.5)
  x <- rowSums(O[row, , drop = FALSE] * L[col, , drop = FALSE]) - a[col] + rnorm(nnz)
  data.frame(row = row, col = col, value = as.numeric(x > 0))
}

sweep_seconds <- function(long, n, p, threads, sweeps = 10) {
  utils::capture.output(r <- mmBPFA_sampler(long, H = 30, K_start = 10, warmup = sweeps, iter = 1, thin = 1,
                                            threads = threads, max_memory_gb = 64))
  r[[1]]$seconds_per_sweep        # timed inside the chain, so data preparation is excluded
}

grid <- expand.grid(np = c("200x100", "2000x1000", "20000x10000"), sparsity = c(0, 0.9, 0.99, 0.999),
                    threads = c(1, 8), stringsAsFactors = FALSE)
grid$n <- as.integer(sub("x.*", "", grid$np)); grid$p <- as.integer(sub(".*x", "", grid$np))
grid <- grid[grid$n * grid$p * (1 - grid$sparsity) <= 2e7 & grid$n * grid$p * (1 - grid$sparsity) >= 1e3, ]

rows <- list()
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  long <- simulate_long(g$n, g$p, g$sparsity)
  secs <- sweep_seconds(long, g$n, g$p, g$threads)
  nnz <- nrow(long)
  rows[[i]] <- data.frame(n = g$n, p = g$p, sparsity = g$sparsity, threads = g$threads, nnz = nnz,
                          ms_per_sweep = 1000 * secs, ns_per_cell_per_column = 1e9 * secs / (nnz * 30))
  print(rows[[i]], row.names = FALSE)
}
out <- do.call(rbind, rows)
args <- commandArgs(TRUE)
if (length(args)) utils::write.csv(out, args[1], row.names = FALSE)
