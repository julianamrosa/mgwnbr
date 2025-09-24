#' Print method for summary.myclass objects
#'
#' @param x An object of class `summary.mgwnbr_model`.
#' @param ... Further arguments passed to print.
#'
#' @return The object (invisibly).
#' @export

## print do summary ##

print.summary.mgwnbr_model <- function(x, ...){
  cat("Model Summary:\n")
  cat("\nGeneral Bandwidth:\n")
  print(x$general_bandwidth)
  if (!is.null(x$bandwidths)){
    cat("\nCovariate Bandwidths:\n")
    print(x$bandwidths)
  }
  cat("\nEffective Number of Parameters and T Test Information:\n")
  print(x$test_info)
  cat("\nFive Number Summary for Parameter Estimates:\n")
  print(x$est_5numb_summary)
  cat("\nFive Number Summary for Standard Errors:\n")
  print(x$se_5numb_summary)
  cat(paste(paste("\nDeviance: ", x$deviance, ",", sep=""),
            "Full Log-Likelihood:", x$full_log_likelihood))
  cat(paste(paste("\nAIC: ", x$aic, ",", sep=""),
            "AICc:", x$aicc))
  invisible(x)
}
