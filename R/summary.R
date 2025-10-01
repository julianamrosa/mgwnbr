#' Summary method for mgwnbr_model objects
#'
#' @param object An object of class `mgwnbr_model`.
#' @param ... Additional arguments passed to other methods.
#'
#' @return An object of class `summary.mgwnbr_model`.
#' @export

## Método summary ##

summary.mgwnbr_model <- function(object, ...){
  #check if object is of class mgwnbr
  if (!inherits(object, "mgwnbr_model")) {
    stop("Object is not of class 'mgwnbr_model'")
  }
  summary_list <- list()
  #general bandwidth
  summary_list$general_bandwidth <- object$general_bandwidth
  #covariate bandwidths
  if (!is.null(object$band)){
    summary_list$bandwidths <- unlist(object$band[nrow(object$band), -ncol(object$band)])
  }
  #enp and t test info
  summary_list$test_info <- rbind(object$ENP, object$alpha_level_5_pct, object$t_critical)
  rownames(summary_list$test_info) <- c("ENP", "alpha_level_5_pct", "t_critical")
  #5 number summary for estimates
  summary_list$estimates <- object$mgwr_param_estimates
  summary_list$est_5numb_summary <- rbind(apply(summary_list$estimates, 2, min),
                                          apply(summary_list$estimates, 2, quantile, probs=0.25),
                                          apply(summary_list$estimates, 2, median),
                                          apply(summary_list$estimates, 2, quantile, probs=0.75),
                                          apply(summary_list$estimates, 2, max))
  rownames(summary_list$est_5numb_summary) <- c('Min', '1Q', 'Median', '3Q',  'Max')
  #5 number summary for standard errors
  summary_list$std_errors <- object$mgwr_se
  summary_list$se_5numb_summary <- rbind(apply(summary_list$std_errors, 2, min),
                                         apply(summary_list$std_errors, 2, quantile, probs=0.25),
                                         apply(summary_list$std_errors, 2, median),
                                         apply(summary_list$std_errors, 2, quantile, probs=0.75),
                                         apply(summary_list$std_errors, 2, max))
  rownames(summary_list$se_5numb_summary) <- c('Min', '1Q', 'Median', '3Q',  'Max')
  summary_list$deviance <- round(as.numeric(object$measures[1]), 4)
  summary_list$full_log_likelihood <- round(as.numeric(object$measures[2]), 4)
  summary_list$aic <- round(as.numeric(object$measures[5]), 4)
  summary_list$aicc <- round(as.numeric(object$measures[6]), 4)
  #list
  # summary_list <- list(general_bandwidth, bandwidths, test_info,
  #                      est_5numb_summary, se_5numb_summary, deviance,
  #                      full_log_likelihood, aic, aicc)
  class(summary_list) <- "summary.mgwnbr_model"
  return(summary_list)
}
