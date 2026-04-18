# ---------------------------------------------------------------------------- #
# S3 methods for class "ccn"
# ---------------------------------------------------------------------------- #

#' @export
print.ccn <- function(x, ...) {
  cat("\nCCN Bunching Correction Estimator\n")
  cat("Method:", x$method, "\n")

  # Human-readable estimand label.
  estimand_label <- switch(
    x$estimand %||% "beta_global",
    "beta_global" = "Global beta (CCN 2023)",
    "beta(0+)"    = "Local beta(0+) (CCT 2025)",
    x$estimand
  )
  cat("Estimand:", estimand_label, "\n")

  cat("Obs:", x$n_obs, " | Censored:", x$n_censored,
      sprintf(" (%.1f%%)\n", 100 * x$n_censored / x$n_obs))

  # --- CCT branch: use bootstrap SE if available ---
  if (isTRUE(x$method == "cct")) {
    cat("\nTreatment effect:\n")
    beta_val <- unname(x$beta[1])
    se_val <- x$boot_se
    ci_lo <- x$boot_ci[1]
    ci_hi <- x$boot_ci[2]
    se_lab <- if (is.na(se_val)) "      NA" else sprintf("%.5f", se_val)
    cat(sprintf("  %-20s %10.5f  (boot SE: %s)\n",
                names(x$beta)[1], beta_val, se_lab))
    if (!is.na(ci_lo) && !is.na(ci_hi)) {
      cat(sprintf("  95%% bootstrap CI: [%.5f, %.5f]  (B_ok = %d)\n",
                  ci_lo, ci_hi, x$boot_B))
    }
    if (!is.null(x$locscale_used)) {
      cat(sprintf("\nLocation-scale regression: %s\n",
                  if (x$locscale_used == "on") "on (full CDK)" else "off (Remark 4.2 short-circuit)"))
    }
    cat(sprintf("pi_hat: %.5f  |  Delta/pi: %.5f\n",
                x$pi_hat %||% NA_real_, x$Delta_over_pi %||% NA_real_))
    cat(sprintf("Z cells: %d  |  alpha* = %.2f\n", x$n_bins, x$alpha_star))
    cat("\n")
    return(invisible(x))
  }

  has_boot <- !is.null(x$boot_se_beta) &&
    isTRUE((x$boot_B %||% 0L) > 0L)

  cat("\nTreatment effect (beta):\n")
  if (has_boot) {
    cat(sprintf("  %-20s %10s  %10s  %10s\n",
                "", "Estimate", "OLS SE", "Boot SE"))
  }

  # handle name mismatch: in the lm the var is "I", but we label with treatment name
  beta_names <- names(x$beta)
  fit_coefs <- coef(x$fit)
  fit_vcov  <- vcov(x$fit)

  for (i in seq_along(x$beta)) {
    nm <- beta_names[i]
    # find matching coef in the fit
    fit_nm <- if (nm == x$treatment) "I" else nm
    se_val <- if (fit_nm %in% names(fit_coefs)) sqrt(fit_vcov[fit_nm, fit_nm]) else NA
    if (has_boot) {
      bse <- x$boot_se_beta[nm]
      bse_str <- if (is.null(bse) || is.na(bse)) "      NA" else sprintf("%10.5f", bse)
      cat(sprintf("  %-20s %10.5f  %10.5f  %s\n",
                  nm, x$beta[i], se_val, bse_str))
    } else {
      cat(sprintf("  %-20s %10.5f  (SE: %.5f)\n", nm, x$beta[i], se_val))
    }
  }

  if (has_boot && !is.null(x$boot_ci_beta)) {
    cat("\n  95% bootstrap CI (beta):\n")
    for (i in seq_along(x$beta)) {
      nm <- beta_names[i]
      lo <- x$boot_ci_beta[nm, 1L]
      hi <- x$boot_ci_beta[nm, 2L]
      if (!is.na(lo) && !is.na(hi)) {
        cat(sprintf("    %-18s [%10.5f, %10.5f]\n", nm, lo, hi))
      }
    }
  }

  if (!is.null(x$delta) && x$method != "naive") {
    cat("\nCorrection term (delta):\n")
    if (has_boot && !is.null(x$boot_se_delta)) {
      cat(sprintf("  %-20s %10s  %10s  %10s\n",
                  "", "Estimate", "OLS SE", "Boot SE"))
    }
    for (i in seq_along(x$delta)) {
      nm <- names(x$delta)[i]
      se_val <- if (nm %in% rownames(fit_vcov)) sqrt(fit_vcov[nm, nm]) else NA
      if (has_boot && !is.null(x$boot_se_delta)) {
        bse <- x$boot_se_delta[nm]
        bse_str <- if (is.null(bse) || is.na(bse)) "      NA" else sprintf("%10.5f", bse)
        cat(sprintf("  %-20s %10.5f  %10.5f  %s\n",
                    nm, x$delta[i], se_val, bse_str))
      } else {
        cat(sprintf("  %-20s %10.5f  (SE: %.5f)\n", nm, x$delta[i], se_val))
      }
    }
    if (has_boot && !is.null(x$boot_ci_delta)) {
      cat("\n  95% bootstrap CI (delta):\n")
      for (i in seq_along(x$delta)) {
        nm <- names(x$delta)[i]
        lo <- x$boot_ci_delta[nm, 1L]
        hi <- x$boot_ci_delta[nm, 2L]
        if (!is.na(lo) && !is.na(hi)) {
          cat(sprintf("    %-18s [%10.5f, %10.5f]\n", nm, lo, hi))
        }
      }
    }
  }

  if (has_boot) {
    cat(sprintf("\nBootstrap: pairs, B_ok = %d\n", x$boot_B))
  }

  if (x$method != "naive") {
    cat(sprintf("\nInstrument bins: %d", x$n_bins))
    if (x$n_switched > 0) {
      cat(sprintf(" | Fallback bins: %d", x$n_switched))
    }
    cat(sprintf("\nMean censored expectation: %.4f\n", x$cens_exp_mean))
  }
  cat("\n")
  invisible(x)
}


#' @export
summary.ccn <- function(object, ...) {
  cat("\nCCN Bunching Correction Estimator\n")
  cat(rep("-", 50), "\n", sep = "")
  cat("Call:\n")
  print(object$call)
  cat("\nMethod:", object$method, "\n")
  cat("Observations:", object$n_obs, "\n")
  cat("Censored:", object$n_censored,
      sprintf("(%.1f%%)\n", 100 * object$n_censored / object$n_obs))

  if (object$method != "naive") {
    cat("Instrument bins:", object$n_bins, "\n")
    if (object$n_switched > 0) {
      cat("Bins with fallback:", object$n_switched, "\n")
    }
    cat(sprintf("Mean censored expectation: %.4f\n", object$cens_exp_mean))
  }

  cat("\n")
  cat("--- Underlying regression ---\n")
  print(summary(object$fit))
  invisible(object)
}


#' @export
coef.ccn <- function(object, ...) {
  if (is.null(object$delta)) return(object$beta)
  c(object$beta, object$delta)
}


#' @export
plot.ccn <- function(x, ...) {
  discontinuity_plot(
    outcome    = x$outcome,
    treatment  = x$treatment,
    data       = x$data,
    instrument = x$instrument,
    ...
  )
}
