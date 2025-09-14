#' Print method for mgwnbr_model objects
#'
#' @param object An object of class `mgwnbr_model`.
#' @param ... Additional arguments passed to other methods.
#'
#' @return The object (invisibly).
#' @export

## Método print ##

print.mgwnbr_model <- function(object, ...) {
  cat("General Bandwidth:\n")
  print(object$general_bandwidth)
  if (!is.null(object$band)){
    cat("\nBandwidths for each Golden Section Search Iteration:\n")
    print(object$band)
  }
  cat("\nFitted Values:\n")
  print(object$fitted_values)
  cat("\nGoodness of Fit Measures:\n")
  print(object$measures)
  cat("\nEffective Number of Parameters:\n")
  print(object$ENP)
  cat("\nParameter Estimates:\n")
  print(object$mgwr_param_estimates)
  cat("\nAlpha Levels for Parameter Significance Tests:\n")
  print(object$alpha_level_5_pct)
  cat("\nCritical Values for Parameter Significance Tests:\n")
  print(object$t_critical)
  cat("\nStandard Errors for Parameter Estimates:\n")
  print(object$mgwr_se)
  cat("\nGlobal Parameter Estimates:\n")
  print(object$global_param_estimates)
  cat("\nDenominator Degrees of Freedom for the T Tests:\n")
  print(object$t_test_dfs)
  cat("\nGoodness of Fit Statistics for the Global Model:\n")
  print(object$global_measures)
  invisible(object)
}
