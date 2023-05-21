#' Kaplan-Meier (KM) estimate at specific time point
#'
#' @param time Character vector. Name of time variable. 
#' @param status Character vector. Name of status variable. 
#' @param t_risk Follow-up value to calculate cases, controls. See details. 
#' 
#' @return KM estimate at specific time point
#'
#' @references Pencina MJ, D'Agostino RB Sr, Steyerberg EW. Extensions of net reclassification 
#'  improvement calculations to measure usefulness of new biomarkers. Stat Med. 2011;30(1):11-21
#' @references Inoue E (2018). nricens: NRI for Risk Prediction Models with Time to Event and Binary 
#'  Response Data. R package version 1.6, <https://CRAN.R-project.org/package=nricens>.
#'
#' @seealso \code{\link{km_fit}}
#' @author Martijn Heymans, 2023
#'
#' @export
km_fit <- function(time, status, t_risk){
  fit_km <-
    survfit(Surv(time, status) ~ 1)
  km_t <-
    stepfun(fit_km$time, c(1, fit_km$surv))(t_risk)
  return(km_t)
}

