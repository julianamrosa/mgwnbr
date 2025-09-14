#' Print method for summary.myclass objects
#'
#' @param object An object of class `summary.mgwnbr_model`.
#' @param ... Further arguments passed to print.
#'
#' @return The object (invisibly).
#' @export

## print do summary ##

print.summary.mgwnbr_model <- function(object, ...){
  cat("Model Summary:\n")
  cat("\nGeneral Bandwidth:\n")
  print(object$general_bandwidth)
  if (!is.null(object$bandwidths)){
    cat("\nCovariate Bandwidths:\n")
    print(object$bandwidths)
  }
  cat("\nEffective Number of Parameters and T Test Information:\n")
  print(object$test_info)
  cat("\nFive Number Summary for Parameter Estimates:\n")
  print(object$est_5numb_summary)
  cat("\nFive Number Summary for Standard Errors:\n")
  print(object$se_5numb_summary)
  cat(paste(paste("\nDeviance: ", object$deviance, ",", sep=""),
            "Full Log-Likelihood:", object$full_log_likelihood))
  cat(paste(paste("\nAIC: ", object$aic, ",", sep=""),
            "AICc:", object$aicc))
  invisible(object)
}
