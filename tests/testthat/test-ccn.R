# ---------------------------------------------------------------------------- #
# Tests for the ccn() estimator
#
# DGP: x* = mu + gamma*z + e, x = max(x*, 0), y = alpha + beta*x* + u
# Naive OLS of y on x is biased because x censors x* at zero.
# The CCN correction should recover beta closer to the truth.
# ---------------------------------------------------------------------------- #

test_that("naive estimator runs and returns ccn object", {
  set.seed(1)
  n <- 2000
  z <- runif(n, 0, 10)
  x_star <- 2 + 0.5 * z + rnorm(n)
  x <- pmax(x_star, 0)
  y <- 1 + 0.8 * x_star + rnorm(n)
  df <- data.frame(y = y, x = x, z = z)

  fit <- ccn("y", "x", "z", df, method = "naive")

  expect_s3_class(fit, "ccn")
  expect_equal(fit$method, "naive")
  expect_true(!is.null(fit$beta))
  expect_true(length(fit$beta) == 1)
  expect_null(fit$delta)
})

test_that("uniform correction runs", {
  set.seed(2)
  n <- 2000
  z <- runif(n, 0, 10)
  x_star <- 2 + 0.5 * z + rnorm(n)
  x <- pmax(x_star, 0)
  y <- 1 + 0.8 * x_star + rnorm(n)
  df <- data.frame(y = y, x = x, z = z)

  fit <- ccn("y", "x", "z", df, method = "uniform")

  expect_s3_class(fit, "ccn")
  expect_equal(fit$method, "uniform")
  expect_true(!is.null(fit$delta))
  expect_true(!is.na(fit$cens_exp_mean))
})

test_that("tobit correction runs", {
  set.seed(3)
  n <- 2000
  z <- runif(n, 0, 10)
  x_star <- 2 + 0.5 * z + rnorm(n)
  x <- pmax(x_star, 0)
  y <- 1 + 0.8 * x_star + rnorm(n)
  df <- data.frame(y = y, x = x, z = z)

  fit <- ccn("y", "x", "z", df, method = "tobit")

  expect_s3_class(fit, "ccn")
  expect_equal(fit$method, "tobit")
})

test_that("het_tobit correction runs", {
  set.seed(30)
  n <- 2000
  z <- runif(n, 0, 10)
  x_star <- 2 + 0.5 * z + rnorm(n)
  x <- pmax(x_star, 0)
  y <- 1 + 0.8 * x_star + rnorm(n)
  df <- data.frame(y = y, x = x, z = z)

  fit <- ccn("y", "x", "z", df, method = "het_tobit")

  expect_s3_class(fit, "ccn")
  expect_equal(fit$method, "het_tobit")
})

test_that("symmetric correction runs", {
  set.seed(4)
  n <- 3000
  z <- runif(n, 0, 10)
  x_star <- 2 + 0.5 * z + rnorm(n)
  x <- pmax(x_star, 0)
  y <- 1 + 0.8 * x_star + rnorm(n)
  df <- data.frame(y = y, x = x, z = z)

  fit <- ccn("y", "x", "z", df, method = "symmetric", n_bins = 10)

  expect_s3_class(fit, "ccn")
  expect_equal(fit$method, "symmetric")
})

test_that("corrected estimate is closer to truth than naive", {
  # DGP: y = 1 + 0.8*x + 0.5*x_star + noise
  # True beta = 0.8 (direct effect of censored x)
  # True delta = 0.5 (effect through latent x_star)
  # Naive OLS on x is biased upward (picks up delta through x=x_star correlation)
  set.seed(10)
  n <- 5000
  z <- runif(n, 0, 10)
  x_star <- -1 + 0.3 * z + rnorm(n)   # substantial censoring at 0 (~30%)
  x <- pmax(x_star, 0)
  y <- 1 + 0.8 * x + 0.5 * x_star + rnorm(n)
  df <- data.frame(y = y, x = x, z = z)

  naive <- ccn("y", "x", "z", df, method = "naive")
  corrected <- ccn("y", "x", "z", df, method = "uniform", n_bins = 10)

  # naive should be biased away from 0.8 (toward 1.3 = 0.8 + 0.5)
  naive_bias <- abs(naive$beta[1] - 0.8)
  corr_bias  <- abs(corrected$beta[1] - 0.8)

  expect_true(corr_bias < naive_bias,
              info = sprintf("Corrected bias (%.4f) should be less than naive (%.4f)",
                             corr_bias, naive_bias))
})

test_that("heterogeneous beta works", {
  set.seed(5)
  n <- 3000
  z <- runif(n, 0, 10)
  w <- rbinom(n, 1, 0.5)
  x_star <- 2 + 0.5 * z + rnorm(n)
  x <- pmax(x_star, 0)
  y <- 1 + 0.8 * x_star + 0.3 * x_star * w + rnorm(n)
  df <- data.frame(y = y, x = x, z = z, w = w)

  fit <- ccn("y", "x", "z", df, method = "uniform", het_beta = "w", n_bins = 10)

  expect_true(length(fit$beta) == 2)
  expect_true("I_w" %in% names(fit$beta))
})

test_that("print method works without error", {
  set.seed(6)
  n <- 2000
  z <- runif(n, 0, 10)
  x_star <- -0.5 + 0.3 * z + rnorm(n)
  x <- pmax(x_star, 0)
  y <- 1 + 0.8 * x_star + rnorm(n)
  df <- data.frame(y = y, x = x, z = z)

  fit <- ccn("y", "x", "z", df, method = "uniform", n_bins = 10)
  expect_output(print(fit), "CCN Bunching Correction Estimator")
})

test_that("coef method returns named vector", {
  set.seed(7)
  n <- 2000
  z <- runif(n, 0, 10)
  x_star <- -0.5 + 0.3 * z + rnorm(n)
  x <- pmax(x_star, 0)
  y <- 1 + 0.8 * x_star + rnorm(n)
  df <- data.frame(y = y, x = x, z = z)

  fit <- ccn("y", "x", "z", df, method = "uniform", n_bins = 10)
  co <- coef(fit)
  expect_true(is.numeric(co))
  expect_true(length(co) >= 2)  # at least beta + delta
})

test_that("missing variable gives informative error", {
  df <- data.frame(y = 1:10, x = 1:10)
  expect_error(ccn("y", "x", "z", df), "not found in data")
})

test_that("n_bins = 1 and constant instrument do not crash", {
  set.seed(25)
  n <- 2000
  z <- runif(n, 0, 10)
  x_star <- -0.5 + 0.3 * z + rnorm(n)
  x <- pmax(x_star, 0)
  y <- 1 + 0.8 * x_star + rnorm(n)
  df <- data.frame(y = y, x = x, z = z)

  # n_bins = 1: every method
  fit1 <- ccn("y", "x", "z", df, method = "naive", n_bins = 1)
  expect_s3_class(fit1, "ccn")

  fit2 <- ccn("y", "x", "z", df, method = "uniform", n_bins = 1)
  expect_s3_class(fit2, "ccn")

  fit3 <- ccn("y", "x", "z", df, method = "tobit", n_bins = 1)
  expect_s3_class(fit3, "ccn")

  fit4 <- ccn("y", "x", "z", df, method = "het_tobit", n_bins = 1)
  expect_s3_class(fit4, "ccn")

  fit5 <- ccn("y", "x", "z", df, method = "symmetric",
              swap = "tobit", n_bins = 1)
  expect_s3_class(fit5, "ccn")

  fit6 <- ccn("y", "x", "z", df, method = "symmetric",
              swap = "het_tobit", n_bins = 1)
  expect_s3_class(fit6, "ccn")

  # constant instrument
  df$z_const <- 5
  fit7 <- ccn("y", "x", "z_const", df, method = "naive")
  expect_s3_class(fit7, "ccn")

  fit8 <- ccn("y", "x", "z_const", df, method = "tobit")
  expect_s3_class(fit8, "ccn")
})

test_that("discontinuity_plot returns ggplot with numeric p-value", {
  set.seed(8)
  n <- 3000
  z <- runif(n, 0, 10)
  x_star <- -0.5 + 0.3 * z + rnorm(n)
  x <- pmax(x_star, 0)
  y <- 1 + 0.8 * x_star + rnorm(n)
  df <- data.frame(y = y, x = x, z = z)

  p <- discontinuity_plot("y", "x", df)
  expect_s3_class(p, "gg")
  # title should contain a numeric p-value, not "NA"
  title_text <- p$labels$title
  expect_false(grepl("NA", title_text),
               info = "p-value should not be NA")
})

test_that("plot.ccn() works on a fitted ccn object", {
  set.seed(20)
  n <- 3000
  z <- runif(n, 0, 10)
  x_star <- -0.5 + 0.3 * z + rnorm(n)
  x <- pmax(x_star, 0)
  y <- 1 + 0.8 * x_star + rnorm(n)
  df <- data.frame(y = y, x = x, z = z)

  fit <- ccn("y", "x", "z", df, method = "symmetric", n_bins = 10)
  p <- plot(fit)
  expect_s3_class(p, "gg")
})

test_that("residualized discontinuity plot uses instrument", {
  set.seed(21)
  n <- 5000
  z <- runif(n, 0, 10)
  x_star <- -0.5 + 0.3 * z + rnorm(n)
  x <- pmax(x_star, 0)
  # w correlated with z so residualizing on z changes the curve
  w <- 0.5 * z + rnorm(n)
  y <- 1 + 0.8 * x_star + 3 * w + rnorm(n)
  df <- data.frame(y = y, x = x, z = z, w = w)

  # without instrument in conditioning set
  p1 <- discontinuity_plot("y", "x", df,
                           controls = "w",
                           residualize = TRUE)
  # with instrument in conditioning set
  p2 <- discontinuity_plot("y", "x", df,
                           controls = "w",
                           instrument = "z",
                           residualize = TRUE)
  expect_s3_class(p1, "gg")
  expect_s3_class(p2, "gg")

  # extract the smooth-curve data from the geom_line layer
  # (layer 2 is geom_line in the plot construction)
  line1 <- ggplot2::ggplot_build(p1)$data[[2]]
  line2 <- ggplot2::ggplot_build(p2)$data[[2]]

  # the y-values of the smooth curve should differ
  # because conditioning on z changes the residualized outcome
  expect_false(
    isTRUE(all.equal(line1$y, line2$y, tolerance = 1e-3)),
    info = "Residualized curves should differ when instrument is added"
  )
})

test_that("swap='tobit' and swap='het_tobit' produce different results", {
  # DGP with heavy censoring so fallback is actually triggered
  set.seed(22)
  n <- 5000
  z <- runif(n, 0, 5)
  # ~60% censoring in some bins to trigger the fallback
  x_star <- -2 + 0.5 * z + rnorm(n)
  x <- pmax(x_star, 0)
  y <- 1 + 0.8 * x_star + rnorm(n)
  df <- data.frame(y = y, x = x, z = z)

  fit_tobit <- ccn("y", "x", "z", df,
                   method = "symmetric", swap = "tobit",
                   n_bins = 5)
  fit_het   <- ccn("y", "x", "z", df,
                   method = "symmetric", swap = "het_tobit",
                   n_bins = 5)

  # if any bins were switched, the two should differ
  if (fit_tobit$n_switched > 0 && fit_het$n_switched > 0) {
    expect_false(
      isTRUE(all.equal(fit_tobit$beta[1], fit_het$beta[1])),
      info = paste0(
        "With switched bins, pooled and per-bin Tobit ",
        "should give different beta estimates"
      )
    )
  }
})

test_that("top-level 'tobit' and 'het_tobit' are distinct estimators", {
  set.seed(31)
  n <- 5000
  z <- runif(n, 0, 5)
  x_star <- -2 + 0.5 * z + rnorm(n)
  x <- pmax(x_star, 0)
  y <- 1 + 0.8 * x_star + rnorm(n)
  df <- data.frame(y = y, x = x, z = z)

  fit_tobit <- ccn("y", "x", "z", df, method = "tobit", n_bins = 5)
  fit_het   <- ccn("y", "x", "z", df, method = "het_tobit", n_bins = 5)

  expect_false(
    isTRUE(all.equal(fit_tobit$beta[1], fit_het$beta[1])),
    info = "Pooled Tobit and per-bin Tobit should not collapse to the same estimator"
  )
})
