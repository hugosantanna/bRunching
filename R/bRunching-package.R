#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats lm predict quantile sd model.matrix pnorm dnorm pt coef
#'   vcov residuals weighted.mean ecdf aggregate as.formula loess median
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar geom_line
#'   geom_ribbon labs theme_minimal scale_x_continuous
#' @importFrom survival survreg Surv
## usethis namespace: end
NULL

# ggplot2 .data pronoun
utils::globalVariables(".data")
