#' Cross-check step 1: generate data and run bRunching R estimates.
#'
#' Writes:
#'   test_data.csv  — seed-fixed simulated data (y, x, z)
#'   r_results.csv  — R estimates from bRunching::ccn() for each method

suppressPackageStartupMessages({
  library(bRunching)
})

set.seed(12345)
n <- 5000
z <- runif(n, 0, 10)
x_star <- -1 + 0.3 * z + rnorm(n)
x <- pmax(x_star, 0)
y <- 1 + 0.8 * x_star + rnorm(n)
df <- data.frame(y = y, x = x, z = z)

out_dir <- getwd()

write.csv(df, file.path(out_dir, "test_data.csv"), row.names = FALSE)

n_bins <- 10L
fit <- function(m) ccn("y", "x", "z", df, method = m, n_bins = n_bins)

r_tab <- rbind(
  data.frame(method = "naive",     beta = fit("naive")$beta[[1]],     delta = NA_real_),
  data.frame(method = "uniform",   beta = fit("uniform")$beta[[1]],   delta = fit("uniform")$delta[[1]]),
  data.frame(method = "het_tobit", beta = fit("het_tobit")$beta[[1]], delta = fit("het_tobit")$delta[[1]]),
  data.frame(method = "symmetric", beta = fit("symmetric")$beta[[1]], delta = fit("symmetric")$delta[[1]])
)

write.csv(r_tab, file.path(out_dir, "r_results.csv"), row.names = FALSE)

cat("R results (n =", n, ", bins =", n_bins, ", true beta = 0.8):\n")
print(r_tab, row.names = FALSE)
