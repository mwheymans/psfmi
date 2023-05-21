#' Risk calculation at specific time point for Cox model
#'
#' @param mod a Cox regression model object.
#' @param t_risk Follow-up value to calculate cases, controls. See details. 
#' 
#' @return Cox regression Risk estimates at specific time point.
#'
#' @references Pencina MJ, D'Agostino RB Sr, Steyerberg EW. Extensions of net reclassification 
#'  improvement calculations to measure usefulness of new biomarkers. Stat Med. 2011;30(1):11-21
#' @references Inoue E (2018). nricens: NRI for Risk Prediction Models with Time to Event and Binary 
#'  Response Data. R package version 1.6, <https://CRAN.R-project.org/package=nricens>.
#'
#' @seealso \code{\link{nri_cox}}
#' 
#' @author Martijn Heymans, 2023
#'
#' @export
risk_coxph <-
  function(mod, t_risk) {
    bh <-
      basehaz(mod)
    h0 <-
      approx(bh$time, bh$hazard, t_risk)$y
    risk <-
      1 - exp( - h0 * exp(mod$linear.predictors ) )
    return(risk)
  }

