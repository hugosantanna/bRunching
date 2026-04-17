#' Monte Carlo: shift bin-1 Y intercept by alpha_1 in {-2, 0, 2}.
#'
#' DGP matches bRunching's own unit test (tests/testthat/test-ccn.R L89-111):
#'   y = 1 + alpha_1 * 1{bin == 1} + 0.8*x + 0.5*x_star + noise
#' True CCN parameters: beta = 0.8, delta = 0.5.
#'
#' Referee prediction: R (saturated bin FE) -> beta_hat ~ 0.8 for all alpha_1.
#' Stata .ado (i.Z_k + noconstant) -> beta_hat biased linearly in alpha_1.

suppressPackageStartupMessages({ library(bRunching) })

run_one <- function(seed, alpha_1, n = 5000L, K = 10L) {
  set.seed(seed)
  z <- runif(n, 0, 10)
  brk <- stats::quantile(z, probs = seq(0, 1, length.out = K + 1L))
  bin <- as.integer(cut(z, breaks = brk, include.lowest = TRUE, labels = FALSE))
  x_star <- -1 + 0.3 * z + rnorm(n)
  x      <- pmax(x_star, 0)
  y      <- 1 + alpha_1 * (bin == 1L) + 0.8 * x + 0.5 * x_star + rnorm(n)
  df <- data.frame(y = y, x = x, z = z)
  meths <- c("naive","uniform","tobit","het_tobit","symmetric")
  sapply(meths, function(m) {
    tryCatch(ccn("y","x","z", df, method = m, n_bins = K)$beta[[1]],
             error = function(e) NA_real_)
  })
}

B      <- 500L
alphas <- c(-2, 0, 2)

out <- do.call(rbind, lapply(seq_along(alphas), function(i) {
  a <- alphas[i]
  cat(sprintf("R Monte Carlo: alpha_1 = %+g  ... ", a))
  t0 <- Sys.time()
  mat <- t(sapply(seq_len(B), function(s) run_one(s + 1000L * i, a)))
  cat(sprintf("done (%.1fs)\n", as.numeric(Sys.time() - t0, units = "secs")))
  data.frame(alpha_1 = a, rep = seq_len(B), mat, check.names = FALSE)
}))

write.csv(out, "mc_r_results.csv", row.names = FALSE)

summarize <- function(x) c(mean = mean(x, na.rm = TRUE),
                            bias = mean(x, na.rm = TRUE) - 0.8,
                            sd   = stats::sd(x, na.rm = TRUE),
                            rmse = sqrt(mean((x - 0.8)^2, na.rm = TRUE)))

cat("\n=== R Monte Carlo summary (true beta = 0.8, B =", B, ") ===\n")
meths <- c("naive","uniform","tobit","het_tobit","symmetric")
for (a in alphas) {
  sub <- out[out$alpha_1 == a, meths]
  tab <- t(sapply(sub, summarize))
  cat(sprintf("\n--- alpha_1 = %+g ---\n", a))
  print(round(tab, 4))
}
