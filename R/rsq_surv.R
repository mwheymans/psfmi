#' R-square calculation for Cox regression models
#'
#' @param fitobj a Cox regression model object of "coxph"
#'
#' @return The value for the explained variance.
#'
#' @references F. Harrell. Regression Modeling Strategies. With Applications to
#'  Linear Models, Logistic and Ordinal Regression, and Survival Analysis. 2nd Edition.
#'  Springer, New York, NY, 2015.
#'
#' @seealso \code{\link{pool_performance}}
#' @author Martijn Heymans, 2021
#'
#' @export
rsq_surv <- function(fitobj){
  n <- fitobj$n
  ltest <- -2 * (fitobj$loglik[1] - fitobj$loglik[2])
  R2_max <- 1 - exp(2 * fitobj$loglik[1]/n)
  rsq <- (1 - exp(-ltest/n))/R2_max
  return(rsq)
}

