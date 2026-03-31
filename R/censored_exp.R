# ---------------------------------------------------------------------------- #
# Censored-expectation estimators
#
# Each function takes treatment vector I, bin assignments z, and a censoring
# indicator, and returns a vector of censored expectations (one per obs).
# ---------------------------------------------------------------------------- #

#' @keywords internal
cens_exp_uniform <- function(I, z, cens) {
  bins <- unique(z)
  out <- rep(NA_real_, length(I))

  for (b in bins) {
    idx <- which(z == b)
    pos <- I[idx] > 0
    if (sum(pos) == 0 || all(pos)) next
    pos_mean  <- mean(I[idx][pos])
    bunching  <- mean(cens[idx])
    out[idx]  <- -pos_mean * (bunching / (1 - bunching))
  }
  out
}


#' @keywords internal
cens_exp_tobit <- function(I, z, cens) {
  bins <- unique(z)
  out <- rep(NA_real_, length(I))
  dropped <- integer(0)

  for (b in bins) {
    idx <- which(z == b)
    I_b <- I[idx]

    # need variation and at least some censored + uncensored obs
    if (length(unique(I_b)) < 2 || all(I_b == 0) || all(I_b > 0)) {
      dropped <- c(dropped, idx)
      next
    }

    fit <- tryCatch(
      survival::survreg(
        survival::Surv(I_b, I_b > 0, type = "left") ~ 1,
        dist = "gaussian"
      ),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      dropped <- c(dropped, idx)
      next
    }

    mu    <- coef(fit)
    sigma <- fit$scale
    # E[I | I <= 0] = mu - sigma * dnorm(mu/sigma) / pnorm(-mu/sigma)
    # (truncated normal expectation below zero)
    ratio <- stats::dnorm(mu / sigma) / stats::pnorm(-mu / sigma)
    e_cens <- mu - sigma * ratio
    out[idx] <- e_cens
  }

  attr(out, "dropped") <- dropped
  out
}


#' @keywords internal
cens_exp_tobit_pooled <- function(I, z, cens) {
  # Pooled Tobit with bin dummies (Stata: tobit I dumX_k_*, ll(0) noconstant)
  n <- length(I)
  out <- rep(NA_real_, n)
  dropped <- integer(0)

  if (all(I == 0) || all(I > 0)) {
    attr(out, "dropped") <- seq_len(n)
    return(out)
  }

  z_f <- factor(z)
  single_bin <- nlevels(z_f) < 2

  fit <- tryCatch({
    if (single_bin) {
      # intercept-only Tobit (no bin dummies)
      survival::survreg(
        survival::Surv(I, I > 0, type = "left") ~ 1,
        dist = "gaussian"
      )
    } else {
      bin_dummies <- stats::model.matrix(~ z_f)[, , drop = FALSE]
      df_tobit <- data.frame(I = I, bin_dummies, check.names = FALSE)
      survival::survreg(
        survival::Surv(I, I > 0, type = "left") ~ . - 1,
        data = df_tobit,
        dist = "gaussian"
      )
    }
  }, error = function(e) NULL)

  if (is.null(fit)) {
    attr(out, "dropped") <- seq_len(n)
    return(out)
  }

  sigma <- fit$scale
  # linear predictor for each obs (mu_i = X_i %*% beta)
  mu <- stats::predict(fit, type = "lp")
  # E[I | I <= 0] = mu - sigma * dnorm(mu/sigma) / pnorm(-mu/sigma)
  ratio <- stats::dnorm(mu / sigma) / stats::pnorm(-mu / sigma)
  e_cens <- mu - sigma * ratio
  out <- e_cens

  attr(out, "dropped") <- dropped
  out
}


#' @keywords internal
cens_exp_symmetric <- function(I, z, cens, swap = "tobit") {
  bins   <- unique(z)
  out    <- rep(NA_real_, length(I))
  switch_bins <- character(0)

  for (b in bins) {
    idx   <- which(z == b)
    I_b   <- I[idx]
    cens_b <- cens[idx]

    # check if >50% censored
    high_cens <- (stats::median(cens_b) == 1)

    if (!high_cens) {
      # --- symmetry-based estimator ---
      hatF_0   <- mean(cens_b)
      op_hatF  <- 1 - hatF_0

      if (op_hatF <= 0 || op_hatF >= 1) {
        switch_bins <- c(switch_bins, as.character(b))
        next
      }

      # empirical CDF of I within this bin
      ecdf_b <- stats::ecdf(I_b)
      sorted <- sort(I_b)

      # find the op_hatF percentile
      # step2: largest I such that ecdf(I) <= op_hatF
      cdf_vals <- ecdf_b(sorted)
      candidates <- sorted[cdf_vals <= op_hatF]
      step2 <- if (length(candidates) > 0) max(candidates) else 0

      # step3: mean of I conditional on I >= step2
      step3 <- mean(I_b[I_b >= step2])

      # step4: censored expectation
      out[idx] <- step2 - step3
    } else {
      switch_bins <- c(switch_bins, as.character(b))
    }
  }

  # fallback for switched bins
  if (length(switch_bins) > 0) {
    if (swap == "tobit") {
      # pooled Tobit with bin dummies (matches Stata: tobit I dumX_k_*, ll(0))
      pooled_exp <- cens_exp_tobit_pooled(I, z, cens)
      for (b in switch_bins) {
        idx <- which(z == as.numeric(b))
        dropped_p <- attr(pooled_exp, "dropped")
        if (length(dropped_p) > 0 && all(idx %in% dropped_p)) next
        out[idx] <- pooled_exp[idx]
      }
    } else if (swap == "het_tobit") {
      # per-bin Tobit, only for switched bins
      het_exp <- cens_exp_tobit(I, z, cens)
      for (b in switch_bins) {
        idx <- which(z == as.numeric(b))
        dropped_h <- attr(het_exp, "dropped")
        if (length(dropped_h) > 0 && all(idx %in% dropped_h)) next
        out[idx] <- het_exp[idx]
      }
    }
  }

  attr(out, "n_switched") <- length(switch_bins)
  out
}
