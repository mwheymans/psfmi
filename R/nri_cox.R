#' Net Reclassification Index for Cox Regression Models
#'
#' \code{nri_cox} Net Reclassification Index for Cox Regression Models
#' 
#' @param data Data frame with relevant predictors
#' @param formula0 A formula object to specify the reference model as normally used by glm.
#'   See under "Details" and "Examples" how these can be specified. 
#' @param formula1 A formula object to specify the new model as normally used by glm.
#' @param t_risk Follow-up value to calculate cases, controls. See details. 
#' @param cutoff A numerical vector that defines the outcome probability cutoff values.
#' @param B A logical  scalar. If TRUE bootstrap confidence intervals are calculated, if FALSE only
#'  the NRI estimates are reported.
#' @param nboot A numerical scalar. Number of bootstrap samples to derive the percentile bootstrap
#'  confidence intervals. Default is 10.
#'
#' @details 
#'  A typical formula object has the form \code{Outcome ~ terms}. Categorical variables has to
#'  be defined as \code{Outcome ~ factor(variable)}, restricted cubic spline variables as
#'  \code{Outcome ~ rcs(variable, 3)}. Interaction terms can be defined as
#'  \code{Outcome ~ variable1*variable2} or \code{Outcome ~ variable1 + variable2 + variable1:variable2}.
#'  All variables in the terms part have to be separated by a "+". If a formula
#'  object is used set predictors, cat.predictors, spline.predictors or int.predictors
#'  at the default value of NULL.
#'  
#'  Follow-up for which cases nd controls are determined. For censored cases before this follow-up 
#'  the expected risk of being a case is calculated by using the Kaplan-Meier value to calculate
#'  the expected number of cases.These expected numbers are used to calculate the NRI proportions 
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
#'   lbpmicox1 <- subset(psfmi::lbpmicox, Impnr==1) # extract one dataset
#'   risk_est <- nri_cox(data=lbpmicox1, formula0 = Surv(Time, Status) ~ Duration + Pain,
#'                    formula1 = Surv(Time, Status) ~ Duration + Pain + Function + Radiation,
#'                    t_risk = 80, cutoff=c(0.45), B=TRUE, nboot=10)
#' 
#' @export
nri_cox <- function(data,
                    formula0,
                    formula1,
                    t_risk,
                    cutoff,
                    B=FALSE,
                    nboot=10){
  
  form <-
    terms(formula0)
  form_vars <-
    attr(form, "term.labels")
  if(is_empty(form_vars))
    stop("\n", "No predictors defined, model is empty")
  time <-
    as.character(attr(form, "variables")[[2]][[2]])
  status <-
    as.character(attr(form, "variables")[[2]][[3]])
  
  # Check data input
  if(!(is.data.frame(data)))
    stop("Data should be a data frame")
  data <- data.frame(as_tibble(data))
  data <- mutate_if(data, is.factor, ~ as.numeric(as.character(.x)))
  if(!all(data[status]==1 | data[status]==0))
    stop("Status should be a 0 - 1 variable")
  
  # set regression formula
  Y <-
    c(paste0("Surv(", time, ",", status, ")~"))
  
  form0 <- formula0
  form1 <- formula1
  
  fit_cox0 <- coxph(form0, data=data, x=TRUE)
  fit_cox1 <- coxph(form1, data=data, x=TRUE)
  
  p0 <- risk_coxph(fit_cox0, t_risk)
  p1 <- risk_coxph(fit_cox1, t_risk)
  
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
  
  tab_orig <- table(cat_orig, cat_new)
  
  df_temp <- data.frame(cat_orig, cat_new, time=data[,time], status=data[,status])
  df_ctr <- df_temp[df_temp$time>t_risk, ]
  tab_ctr <- table(df_ctr$cat_orig, df_ctr$cat_new)
  
  df_case <- df_temp[df_temp$time<=t_risk & df_temp$status==1, ]
  n_excl <- nrow(data) - nrow(df_ctr) - nrow(df_case)
  tab_case <- table(df_case$cat_orig, df_case$cat_new)
  
  n_total <- nrow(data)
  n_control <- nrow(df_ctr)
  n_cases <- nrow(df_case)
  n_totals <- matrix(c(n_total,
                       n_cases,
                       n_control,
                       n_excl), 4, 1)
  dimnames(n_totals) <- list(c("Total",
                               "Controls",
                               "Cases",
                               "Censored (excluded)"),
                             c("N"))
  
  res_nri <- nri_est(data,
                     p0=p0,
                     p1=p1,
                     time = time,
                     status = status,
                     cutoff = cutoff,
                     t_risk=t_risk)
  res_nri <- data.frame(res_nri)
  names(res_nri) <- "Estimates"
  rownames(res_nri) <- c("Pr(Up - Case)",
                         "Pr(Down - Case)",
                         "Pr(Up - Control)",
                         "Pr(Down - Control)",
                         "NRI Case",
                         "NRI Control",
                         "NRI")
  
  if(B==TRUE){
    boot_data <- bootstraps(data, times = nboot)
    boot_nri <- lapply(boot_data$splits, function(x) {
      x <- as.data.frame(x)
      fit_cox0 <- coxph(form0, data=x, x=TRUE)
      fit_cox1 <- coxph(form1, data=x, x=TRUE)
      
      p0 <- risk_coxph(fit_cox0, t_risk)
      p1 <- risk_coxph(fit_cox1, t_risk)
      res <- nri_est(x, p0, p1, time = time,
                     status = status,
                     cutoff = cutoff,
                     t_risk=t_risk)
    })
    boot_res <- do.call("rbind", boot_nri)
    
    boot_ci <- t(apply(boot_res, 2,
                       function(x)
                         quantile(x, probs=c(0.025, 0.975))))
    res_nri <- cbind(res_nri, boot_ci)
  }
  
  res_km <- km_estimates(data,
                         p0=p0,
                         p1=p1,
                         time = time,
                         status = status,
                         cutoff=cutoff,
                         t_risk=t_risk)
  km_est <- data.frame(matrix(c(res_km), 2, 3, byrow = TRUE))
  names(km_est) <- c("N Up", "N Down", "Total N")
  rownames(km_est) <- c("Cases", "Controls")
  
  res <- list(data=data,
              prob_orig=p0,
              prob_new=p1,
              time=time,
              status=status,
              cutoff=cutoff,
              t_risk=t_risk,
              reclass_totals=tab_orig,
              reclass_cases=tab_case,
              reclass_controls=tab_ctr,
              totals=n_totals,
              km_est=km_est,
              nri_est=res_nri)
  return(res)
}