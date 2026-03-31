#' Discontinuity Plot at the Bunching Point
#'
#' Creates a plot showing the outcome mean at the bunching point (zero) versus
#' a local-linear extrapolation from positive treatment values, with confidence
#' intervals and a p-value for the discontinuity test.
#'
#' @param outcome Character string. Name of the outcome variable in `data`.
#' @param treatment Character string. Name of the treatment variable in `data`.
#' @param data A data.frame.
#' @param controls Character vector of control variable names for
#'   residualization. If `NULL`, uses raw outcome.
#' @param instrument Character string. Name of the instrument for
#'   residualization (used in the conditioning set). Optional.
#' @param bandwidth Bandwidth for local linear regression. If `NULL`, uses
#'   the default from [stats::loess()].
#' @param max_treatment Numeric. Maximum treatment value to include (default 8,
#'   matching the Stata code).
#' @param bin_width Numeric. Width of bins for the scatter plot (default 0.1).
#' @param residualize Logical. If `TRUE` and `controls` is provided,
#'   residualize the outcome on controls before plotting.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' n <- 5000
#' z <- runif(n, 0, 10)
#' x_star <- 2 + 0.5 * z + rnorm(n)
#' x <- pmax(x_star, 0)
#' y <- 1 + 0.8 * x_star + rnorm(n)
#' df <- data.frame(y = y, x = x, z = z)
#' discontinuity_plot("y", "x", df)
#' }
#'
#' @export
discontinuity_plot <- function(outcome,
                               treatment,
                               data,
                               controls = NULL,
                               instrument = NULL,
                               bandwidth = NULL,
                               max_treatment = 8,
                               bin_width = 0.1,
                               residualize = FALSE) {

  S <- data[[outcome]]
  I <- data[[treatment]]

  # filter to max treatment
  keep <- I <= max_treatment & !is.na(S) & !is.na(I)
  S <- S[keep]
  I <- I[keep]
  sub_data <- data[keep, , drop = FALSE]

  # residualize if requested (Stata: reg LHS RHS controls i.Z_K if RHS > 0)
  if (residualize && !is.null(controls)) {
    rhs_terms <- c(treatment, controls)
    if (!is.null(instrument) && instrument %in% names(sub_data)) {
      sub_data[[".z_factor"]] <- factor(sub_data[[instrument]])
      rhs_terms <- c(rhs_terms, ".z_factor")
    }
    ctrl_fmla <- stats::as.formula(
      paste(outcome, "~", paste(rhs_terms, collapse = " + "))
    )
    pos_mask <- sub_data[[treatment]] > 0
    res_fit <- stats::lm(ctrl_fmla, data = sub_data, subset = pos_mask)
    S_res <- stats::residuals(res_fit) +
      coef(res_fit)[treatment] * I[pos_mask]
    S[pos_mask] <- S_res
  }

  # bin the treatment
  breaks <- c(0, seq(bin_width / 1000, max_treatment + bin_width, by = bin_width))
  x_bin <- cut(I, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  bin_mids <- c(0, seq(bin_width, max_treatment, by = bin_width))
  if (length(bin_mids) < max(x_bin, na.rm = TRUE)) {
    bin_mids <- c(bin_mids, rep(NA, max(x_bin, na.rm = TRUE) - length(bin_mids)))
  }

  # compute bin-level statistics
  bin_stats <- stats::aggregate(
    S,
    by = list(bin = x_bin),
    FUN = function(y) c(mean = mean(y, na.rm = TRUE),
                        sd = stats::sd(y, na.rm = TRUE),
                        n = sum(!is.na(y)))
  )
  bin_stats <- data.frame(
    bin  = bin_stats$bin,
    mean = bin_stats$x[, "mean"],
    sd   = bin_stats$x[, "sd"],
    n    = bin_stats$x[, "n"]
  )
  bin_stats$x   <- bin_mids[bin_stats$bin]
  bin_stats$ci_ub <- bin_stats$mean + 1.96 * bin_stats$sd / sqrt(bin_stats$n)
  bin_stats$ci_lb <- bin_stats$mean - 1.96 * bin_stats$sd / sqrt(bin_stats$n)

  # bunching point stats (bin 1, which is x == 0)
  bp <- bin_stats[bin_stats$x == 0, ]

  # --- extrapolation to zero via local linear regression ---
  # Kernel-weighted local linear, matching Stata's lpoly degree(1).
  pos_idx <- I > 0
  pos_df  <- data.frame(yy = S[pos_idx], xx = I[pos_idx])

  bw <- if (!is.null(bandwidth)) bandwidth else stats::quantile(pos_df$xx, 0.25)

  pred_zero <- lpoly_fit(pos_df$xx, pos_df$yy, eval_at = 0, h = bw)
  lpoly_mean <- pred_zero$fit
  lpoly_se   <- pred_zero$se

  # discontinuity test
  actual_mean <- bp$mean[1]
  actual_se   <- bp$sd[1] / sqrt(bp$n[1])
  diff_mean   <- lpoly_mean - actual_mean
  diff_se     <- sqrt(lpoly_se^2 + actual_se^2)
  diff_t      <- abs(diff_mean / diff_se)
  df_val      <- bp$n[1] - 1
  pvalue      <- 2 * stats::pt(-diff_t, df = max(df_val, 1))
  if (is.na(pvalue)) {
    pvalue_lab <- "NA"
  } else if (pvalue < 0.001) {
    pvalue_lab <- "< 0.001"
  } else {
    pvalue_lab <- sprintf("%.3f", pvalue)
  }

  # --- smooth curve via kernel-weighted local linear (matches lpoly) ---
  grid <- seq(0, max_treatment, length.out = 200)
  lpoly_res <- lpoly_fit(pos_df$xx, pos_df$yy, grid, bw)

  smooth_df <- data.frame(
    x     = grid,
    y     = lpoly_res$fit,
    upper = lpoly_res$fit + 1.96 * lpoly_res$se,
    lower = lpoly_res$fit - 1.96 * lpoly_res$se
  )
  smooth_df <- smooth_df[!is.na(smooth_df$y), ]

  # --- plot ---
  p <- ggplot2::ggplot() +
    # CI ribbon for local linear
    ggplot2::geom_ribbon(
      data = smooth_df,
      ggplot2::aes(x = .data$x, ymin = .data$lower, ymax = .data$upper),
      alpha = 0.25, fill = "grey50"
    ) +
    # local linear line
    ggplot2::geom_line(
      data = smooth_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      color = "black"
    ) +
    # bunching point with CI
    ggplot2::geom_point(
      data = bp,
      ggplot2::aes(x = .data$x, y = .data$mean),
      shape = 1, size = 3
    ) +
    ggplot2::geom_errorbar(
      data = bp,
      ggplot2::aes(x = .data$x, ymin = .data$ci_lb, ymax = .data$ci_ub),
      width = 0.1
    ) +
    ggplot2::scale_x_continuous(breaks = 0:max_treatment) +
    ggplot2::labs(
      x = treatment,
      y = outcome,
      title = paste0("P-value of Discontinuity: ", pvalue_lab)
    ) +
    ggplot2::theme_minimal()

  p
}
