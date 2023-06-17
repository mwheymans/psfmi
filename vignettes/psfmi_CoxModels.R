## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------

  library(psfmi)
  pool_coxr <- psfmi_coxr(data=lbpmicox, nimp=5, impvar="Impnr", 
                formula = Surv(Time, Status) ~ Duration + Radiation + Onset + 
                Function + Age + Previous + Tampascale + JobControl + 
                JobDemand + Social + factor(Expect_cat), p.crit=1,
                method="D1", direction = "BW")
  
  pool_coxr$RR_model
  pool_coxr$multiparm
  

## ---- eval=TRUE---------------------------------------------------------------

  library(psfmi)
  pool_coxr <- psfmi_coxr(data=lbpmicox, nimp=5, impvar="Impnr", 
                formula = Surv(Time, Status) ~ Duration + Radiation + Onset + 
                Function + Age + Previous + Tampascale + JobControl + 
                JobDemand + Social + factor(Expect_cat), p.crit=0.05,
                method="D1", direction = "FW")
  
  pool_coxr$RR_model_final
  pool_coxr$multiparm_final
  pool_coxr$predictors_in
  

## -----------------------------------------------------------------------------

  library(psfmi)
  
  pool_coxr <- psfmi_coxr(data=lbpmicox, nimp=5, impvar="Impnr", 
                formula = Surv(Time, Status) ~ Duration + Radiation + Onset + 
                Function + Age + Previous + Tampascale + factor(Expect_cat) +
                factor(Satisfaction) + Tampascale:Radiation + 
                factor(Expect_cat):Tampascale, keep.predictors = "Tampascale",
                p.crit=0.05, method="D1", direction = "FW")
  
  pool_coxr$RR_model_final
  pool_coxr$multiparm_final
  pool_coxr$predictors_in


## -----------------------------------------------------------------------------

  library(psfmi)
  
  pool_coxr <- psfmi_coxr(data=lbpmicox, nimp=5, impvar="Impnr", 
                formula = Surv(Time, Status) ~ Duration + Radiation + Onset + 
                Function + Previous + rcs(Tampascale, 3) + 
                factor(Satisfaction) + rcs(Tampascale, 3):Radiation,  
                keep.predictors = "Tampascale",
                p.crit=0.05, method="D1", direction = "BW")
  
  pool_coxr$RR_model_final
  pool_coxr$multiparm_final
  pool_coxr$predictors_in


## -----------------------------------------------------------------------------

  library(psfmi)
  pool_coxr <- psfmi_coxr(data=lbpmicox, nimp=5, impvar="Impnr", 
                formula = Surv(Time, Status) ~ Duration + Radiation + Onset + 
                Function + Previous + rcs(Tampascale, 3) + 
                factor(Satisfaction) + rcs(Tampascale, 3):Radiation,  
                keep.predictors = "Tampascale",
                p.crit=0.05, method="MPR", direction = "FW")
  
  pool_coxr$RR_model_final
  pool_coxr$multiparm_final
  pool_coxr$predictors_in
  

## -----------------------------------------------------------------------------

  library(psfmi)
  pool_coxr <- psfmi_coxr(data=lbpmicox, nimp=5, impvar="Impnr", 
                formula = Surv(Time, Status) ~ Duration + Onset + 
                Function + Previous + rcs(Tampascale, 3) + 
                factor(Satisfaction) + strata(Radiation), 
                p.crit=0.05, method="MPR", direction = "BW")
  
  pool_coxr$RR_model_final
  pool_coxr$multiparm_final
  pool_coxr$formula_step
  

