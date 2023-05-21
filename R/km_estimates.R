#' Kaplan-Meier risk estimates for Net Reclassification Index analysis
#'
#' \code{km_estimates} Kaplan-Meier risk estimates for Net Reclassification Index analysis 
#'  for Cox Regression Models
#' 
#' @param data Data frame with relevant predictors
#' @param p0 risk outcome probabilities for reference model. 
#' @param p1 risk outcome probabilities for new model.
#' @param time Character vector. Name of time variable. 
#' @param status Character vector. Name of status variable. 
#' @param t_risk Follow-up value to calculate cases, controls. See details. 
#' @param cutoff A numerical vector that defines the outcome probability cutoff values.
#'
#' @details 
#'  Follow-up for which cases nd controls are determined. For censored cases before this follow-up 
#'  the expected risk of being a case is calculated by using the Kaplan-Meier value to calculate
#'  the expected number of cases. These expected numbers are used to calculate the NRI proportions 
#'  but are not shown by function \code{nricens}.
#'   
#'@return An object from which the following objects can be extracted: 
#'  \itemize{
#'  \item  \code{data} dataset. 
#'  \item  \code{prob_orig} outcome risk probabilities at t_risk for reference model.
#'  \item  \code{prob_new} outcome risk probabilities at t_risk for new model.
#'  \item  \code{time} name of time variable.
#'  \item  \code{status} name of status variable.
#'  \item  \code{cutoff} cutoff value for survival probability. 
#'  \item  \code{t_risk} follow-up time used to calculate outcome (risk) probabilities.
#'  \item  \code{reclass_totals} table with total reclassification numbers.
#'  \item  \code{reclass_cases} table with reclassification numbers for cases.
#'  \item  \code{reclass_controls} table with reclassification numbers for controls.
#'  \item  \code{totals} totals of controls, cases, censored cases.
#'  \item  \code{km_est} totals of cases calculated using Kaplan-Meiers risk estimates.
#'  \item  \code{nri_est} reclassification measures.
#' }
#'
#' @references Cook NR, Ridker PM. Advances in measuring the effect of individual predictors of 
#'  cardiovascular risk: the role of reclassification measures. Ann Intern Med. 2009;150(11):795-802.
#' @references Steyerberg EW, Pencina MJ. Reclassification calculations for persons with incomplete 
#'  follow-up. Ann Intern Med. 2010;152(3):195-6 (author reply 196-7).
#' @references Pencina MJ, D'Agostino RB Sr, Steyerberg EW. Extensions of net reclassification 
#'  improvement calculations to measure usefulness of new biomarkers. Stat Med. 2011;30(1):11-21
#' @references Inoue E (2018). nricens: NRI for Risk Prediction Models with Time to Event and Binary 
#'  Response Data. R package version 1.6, <https://CRAN.R-project.org/package=nricens>.
#'  
#' @author Martijn Heymans, 2023
#'
#' @examples
#'   library(survival)
#'   lbpmicox1 <- subset(psfmi::lbpmicox, Impnr==1) # extract dataset
#'   
#'   fit_cox0 <- 
#'       coxph(Surv(Time, Status) ~ Duration + Pain, data=lbpmicox1, x=TRUE)
#'   fit_cox1 <- 
#'       coxph(Surv(Time, Status) ~ Duration + Pain + Function + Radiation, 
#'       data=lbpmicox1, x=TRUE)
#'
#'   p0 <- risk_coxph(fit_cox0, t_risk=80)
#'   p1 <- risk_coxph(fit_cox1, t_risk=80)
#'   
#'   res_km <- km_estimates(data=lbpmicox1,
#'                       p0=p0,
#'                       p1=p1,
#'                       time = "Time",
#'                       status = "Status",
#'                       cutoff=0.45,
#'                       t_risk=80)
#' 
#' @export
km_estimates <- function(data, p0, p1, time, status, t_risk, cutoff){
  
  nctr_km <-
    nrow(data)*km_fit(data[,time], data[,status], t_risk)
  ncas_km <-
    nrow(data)* (1-km_fit(data[,time], data[,status], t_risk))
  
  cat_orig <- cut(p0,
                  breaks = c(0, cutoff, 1),
                  include.lowest = TRUE,
                  right=FALSE,
  )
  cat_new <- cut(p1,
                 breaks = c(0, cutoff, 1),
                 include.lowest = TRUE,
                 right=FALSE,
  )
  
  ####################################################
  
  cat0 <-
    as.numeric(cat_orig)
  cat1 <-
    as.numeric(cat_new)
  
  n_down <-
    cat1-cat0 < 0
  n_up <-
    cat1-cat0 > 0
  
  if(sum(n_up)==0){
    n_up <- n_up_case <- n_up_ctr <- 0
  } else {
    # KM probability
    # Up
    km_up_case <-
      1-km_fit(data[,time][n_up], data[,status][n_up], t_risk)
    km_up_ctr <-
      km_fit(data[,time][n_up], data[,status][n_up], t_risk)
    
    n_up_case <-
      sum(n_up)*km_up_case
    n_up_ctr <-
      sum(n_up)*km_up_ctr
  }
  if(sum(n_down)==0){
    n_down <- n_down_case <- n_down_ctr <- 0
  } else {
    # down
    km_down_case <-
      1-km_fit(data[,time][n_down], data[,status][n_down], t_risk)
    km_down_ctr <-
      km_fit(data[,time][n_down], data[,status][n_down], t_risk)
    
    n_down_case <-
      sum(n_down)*km_down_case
    n_down_ctr <-
      sum(n_down)*km_down_ctr
  }
  
  ###################################################
  
  km_obj <- c(n_up_case,
              n_down_case,
              ncas_km,
              n_up_ctr,
              n_down_ctr,
              nctr_km)
  return(km_obj)
}