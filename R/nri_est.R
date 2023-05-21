#' Calculation of Net Reclassification Index measures
#'
#' \code{nri_est} Calculation of proportion of Reclassified persons and NRI for Cox 
#'  Regression Models
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
#'  \item  \code{prop_up_case} proportion of cases reclassified upwards. 
#'  \item  \code{prop_down_case} proportion of cases reclassified downwards.
#'  \item  \code{prop_up_ctr} proportion of controls reclassified upwards.
#'  \item  \code{prop_down_ctr} proportion of controls reclassified downwards.
#'  \item  \code{nri_plus} proportion reclassified for events.
#'  \item  \code{nri_min} proportion reclassified for nonevents. 
#'  \item  \code{nri} net reclassification improvement.
#' }
#'
#' @references Cook NR, Ridker PM. Advances in measuring the effect of individual predictors of 
#'  cardiovascular risk: the role of reclassification measures. Ann Intern Med. 2009;150(11):795-802.
#' @references Steyerberg EW, Pencina MJ. Reclassification calculations for persons with incomplete 
#'  follow-up. Ann Intern Med. 2010;152(3):195-6; author reply 196-7.
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
#'     coxph(Surv(Time, Status) ~ Duration + Pain, data=lbpmicox1, x=TRUE)
#'   fit_cox1 <- 
#'     coxph(Surv(Time, Status) ~ Duration + Pain + Function + Radiation, 
#'     data=lbpmicox1, x=TRUE)
#'
#'   p0 <- risk_coxph(fit_cox0, t_risk=80)
#'   p1 <- risk_coxph(fit_cox1, t_risk=80)
#'   
#'   nri <- nri_est(data=lbpmicox1,
#'                       p0=p0,
#'                       p1=p1,
#'                       time = "Time",
#'                       status = "Status",
#'                       cutoff=0.45,
#'                       t_risk=80)
#' 
#' @export
nri_est <- function(data, p0, p1, time, status, t_risk, cutoff){
  
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
    n_up <- prop_up_case <- prop_up_ctr <- 0
  } else {
    # KM probability
    # Up
    km_up_case <-
      1-km_fit(data[,time][n_up], data[,status][n_up], t_risk)
    km_up_ctr <-
      km_fit(data[,time][n_up], data[,status][n_up], t_risk)
    
    n_up_case <-
      sum(n_up)*km_up_case
    prop_up_case <-
      n_up_case / ncas_km
    
    n_up_ctr <-
      sum(n_up)*km_up_ctr
    prop_up_ctr <-
      n_up_ctr / nctr_km
  }
  if(sum(n_down)==0){
    n_down <- prop_down_case <- prop_down_ctr <- 0
  } else {
    # down
    km_down_case <-
      1-km_fit(data[,time][n_down], data[,status][n_down], t_risk)
    km_down_ctr <-
      km_fit(data[,time][n_down], data[,status][n_down], t_risk)
    
    n_down_case <-
      sum(n_down)*km_down_case
    prop_down_case <-
      n_down_case / ncas_km
    
    n_down_ctr <-
      sum(n_down)*km_down_ctr
    prop_down_ctr <-
      n_down_ctr / nctr_km
  }
  
  ###################################################
  
  nri_plus <-
    prop_up_case - prop_down_case
  nri_min <-
    prop_down_ctr - prop_up_ctr
  nri <-
    nri_plus + nri_min
  
  nri_obj <- c(prop_up_case,
               prop_down_case,
               prop_up_ctr,
               prop_down_ctr,
               nri_plus,
               nri_min,
               nri)
  return(nri_obj)
}