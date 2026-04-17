# ---------------------------------------------------------------------------- #
# Tests for the CCT-2025 branch of ccn() and its helpers.
#
# DGP (strategy memo, Section 5):
#   Z           ~ sample(1:K, n, replace = TRUE)         # K = 20
#   mu(z)       = -0.5 + 0.15 * z
#   eta         ~ Normal(0, 1)
#   X_star      = mu(Z) + eta
#   X           = pmax(X_star, 0)
#   Y           = 2 + 0.8 * X_star + 1.5 * X * (X > 0) + rnorm(n)
# True beta(0+) = 1.5.
# ---------------------------------------------------------------------------- #

test_that("cct_m_hat recovers m(Z) on the homoskedastic normal DGP", {
  set.seed(1)
  n <- 5000
  K <- 20
  Z <- sample.int(K, n, replace = TRUE)
  mu_z <- -0.5 + 0.15 * Z
  eta  <- rnorm(n)
  x_star <- mu_z + eta
  x <- pmax(x_star, 0)

  sdX <- sd(x)
  m_res <- bRunching:::cct_m_hat(
    X          = x,
    Z          = Z,
    alpha_star = 0.5,
    alpha_grid = seq(0.55, 0.95, by = 0.05),
    alpha_s    = 0.95,
    zeta_0     = 1e-3 * sdX,
    zeta_1     = 1e-3 * sdX,
    locscale   = "auto"
  )

  # True m(z) is mu_z = -0.5 + 0.15 * z.
  z_levels <- sort(unique(Z))
  true_m <- -0.5 + 0.15 * z_levels
  m_hat_z <- m_res$m_hat_by_z[as.character(z_levels)]

  # Tolerance: within 3 * MC SE per cell. SE ~ sd(X) / sqrt(n/K) ~ 1 / sqrt(250) ~ 0.063
  # so 3 * SE ~ 0.19 per cell; allow 0.30 to account for finite-sample
  # slack in the CDK weighted OLS step.
  dev <- abs(m_hat_z - true_m)
  expect_true(max(dev, na.rm = TRUE) < 0.30,
              info = sprintf("max |m_hat - true_m| = %.4f", max(dev, na.rm = TRUE)))
  # mean deviation should be small
  expect_true(mean(dev, na.rm = TRUE) < 0.15,
              info = sprintf("mean |m_hat - true_m| = %.4f", mean(dev, na.rm = TRUE)))
})


test_that("ccn(method = 'cct') recovers beta(0+) = 1.5 on the MC DGP", {
  set.seed(1)
  n <- 5000
  K <- 20
  Z <- sample.int(K, n, replace = TRUE)
  mu_z <- -0.5 + 0.15 * Z
  eta  <- rnorm(n)
  x_star <- mu_z + eta
  x <- pmax(x_star, 0)
  y <- 2 + 0.8 * x_star + 1.5 * x * (x > 0) + rnorm(n)
  df <- data.frame(y = y, x = x, z = Z)

  fit <- ccn("y", "x", "z", df,
             method = "cct",
             boot_B = 100L,
             seed   = 42L)

  expect_s3_class(fit, "ccn")
  expect_equal(fit$method, "cct")
  expect_equal(fit$estimand, "beta(0+)")
  expect_equal(names(fit$beta), "beta(0+)")
  expect_true(is.null(fit$delta))
  expect_true(length(fit$beta) == 1L)

  # Tolerance: 3 * bootstrap SE around the truth 1.5
  se <- fit$boot_se
  expect_true(is.finite(se) && se > 0,
              info = sprintf("bootstrap SE = %s", format(se)))
  expect_true(abs(fit$beta - 1.5) < 3 * se,
              info = sprintf("beta_hat = %.4f, SE = %.4f, truth = 1.5",
                             fit$beta, se))

  # bootstrap count: most reps should succeed
  expect_true(fit$boot_B >= 80L,
              info = sprintf("only %d of 100 bootstrap reps succeeded", fit$boot_B))
})


test_that("locscale = 'off' short-circuits to m_hat(z) = median(X | Z = z)", {
  set.seed(99)
  n <- 3000
  K <- 10
  # Shift mu_z up so that every within-z median is strictly positive.
  Z <- sample.int(K, n, replace = TRUE)
  mu_z <- 1 + 0.2 * Z            # all > 0
  eta  <- rnorm(n)
  x_star <- mu_z + eta
  x <- pmax(x_star, 0)           # almost no censoring here

  sdX <- sd(x)
  m_res <- bRunching:::cct_m_hat(
    X          = x,
    Z          = Z,
    alpha_star = 0.5,
    alpha_grid = seq(0.55, 0.95, by = 0.05),
    alpha_s    = 0.95,
    zeta_0     = 1e-3 * sdX,
    zeta_1     = 1e-3 * sdX,
    locscale   = "off"
  )

  # Expected: m_hat(z) = sample median of X within z-cell (type = 1).
  z_levels <- sort(unique(Z))
  medians <- sapply(z_levels, function(zv) stats::quantile(x[Z == zv], 0.5, type = 1, names = FALSE))
  expect_equal(unname(m_res$m_hat_by_z), unname(medians), tolerance = 1e-12)
  expect_equal(m_res$locscale_used, "off")
})
