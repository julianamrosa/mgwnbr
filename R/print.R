#' Print method for mgwnbr_model objects
#'
#' @param x An object of class `mgwnbr_model`.
#' @param ... Additional arguments passed to other methods.
#'
#' @return The object (invisibly).
#' @export

## Método print ##

print.mgwnbr_model <- function(x, ...) {
  cat("General Bandwidth:\n")
  print(x$general_bandwidth)
  if (!is.null(x$band)){
    cat("\nBandwidths for each Golden Section Search Iteration:\n")
    print(x$band)
  }
  cat("\nFitted Values:\n")
  print(x$fitted_values)
  cat("\nGoodness of Fit Measures:\n")
  print(x$measures)
  cat("\nEffective Number of Parameters:\n")
  print(x$ENP)
  cat("\nParameter Estimates:\n")
  print(x$mgwr_param_estimates)
  cat("\nAlpha Levels for Parameter Significance Tests:\n")
  print(x$alpha_level_5_pct)
  cat("\nCritical Values for Parameter Significance Tests:\n")
  print(x$t_critical)
  cat("\nStandard Errors for Parameter Estimates:\n")
  print(x$mgwr_se)
  cat("\nGlobal Parameter Estimates:\n")
  print(x$global_param_estimates)
  cat("\nDenominator Degrees of Freedom for the T Tests:\n")
  print(x$t_test_dfs)
  cat("\nGoodness of Fit Statistics for the Global Model:\n")
  print(x$global_measures)
  invisible(x)
}
