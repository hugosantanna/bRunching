#' Diagnostic: recompute per-bin stats in R and compare to Stata.

df <- read.csv("test_data.csv")
n <- nrow(df); K <- 10L

# match Stata's xtile: equal-frequency bins via quantile breaks
brk <- stats::quantile(df$z, probs = seq(0, 1, length.out = K + 1L))
z_bins <- as.integer(cut(df$z, breaks = brk, include.lowest = TRUE, labels = FALSE))

cens <- as.integer(df$x == 0)
pos  <- df$x > 0

r_stats <- do.call(rbind, lapply(seq_len(K), function(b) {
  idx <- which(z_bins == b)
  pos_idx <- idx[pos[idx]]
  pos_mean <- mean(df$x[pos_idx])
  bunching <- mean(cens[idx])
  cens_exp <- -pos_mean * (bunching / (1 - bunching))
  data.frame(bin = b, pos_mean = pos_mean, bunching = bunching, cens_exp = cens_exp)
}))

cat("=== R per-bin stats ===\n")
print(r_stats, row.names = FALSE, digits = 8)

st <- read.csv("stata_bin_stats.csv")
st <- st[order(st$X_10), ]

cmp <- data.frame(
  bin            = r_stats$bin,
  r_cens_exp     = r_stats$cens_exp,
  stata_cens_exp = st$cens_exp,
  diff           = r_stats$cens_exp - st$cens_exp
)
cat("\n=== R vs Stata cens_exp per bin ===\n")
print(cmp, row.names = FALSE, digits = 8)
cat("\nmax abs diff:", max(abs(cmp$diff)), "\n")

# Now do the actual regression with these in R and compare
reg_term <- r_stats$cens_exp[z_bins] * cens + df$x

bin_f <- factor(z_bins)
bd <- stats::model.matrix(~ 0 + bin_f)
colnames(bd) <- paste0("bin_", seq_len(ncol(bd)))

dfr <- data.frame(S = df$y, I = df$x, bd, reg_term = reg_term,
                  check.names = FALSE)
fit <- stats::lm(S ~ 0 + ., data = dfr)

cat("\n=== R manual regression (match Stata het_uniform) ===\n")
cat(sprintf("beta     = %.6f\n", coef(fit)["I"]))
cat(sprintf("reg_term = %.6f\n", coef(fit)["reg_term"]))

library(bRunching)
pkg_fit <- ccn("y", "x", "z", df, method = "uniform", n_bins = 10)
cat("\n=== bRunching::ccn() uniform ===\n")
cat(sprintf("beta     = %.6f\n", pkg_fit$beta[[1]]))
cat(sprintf("reg_term = %.6f\n", pkg_fit$delta[[1]]))
