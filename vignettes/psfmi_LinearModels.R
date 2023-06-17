## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------

  library(psfmi)
  pool_lm <- psfmi_lm(data=lbpmilr, nimp=5, impvar="Impnr", 
                      formula = Pain ~ Gender + Smoking + 
                      Function + JobControl + JobDemands + SocialSupport, 
                      method="D1")
  
  pool_lm$RR_model
 

## -----------------------------------------------------------------------------

  library(psfmi)
  pool_lm <- psfmi_lm(data=lbpmilr, nimp=5, impvar="Impnr", 
                      formula = Pain ~ Gender + Smoking + 
                      Function + JobControl + JobDemands + SocialSupport, 
                      keep.predictors = "Smoking", method="D1", p.crit=0.05, 
                      direction="BW")
  
  pool_lm$RR_model_final
  pool_lm$multiparm_final
  pool_lm$predictors_out
  

## -----------------------------------------------------------------------------

  library(psfmi)
  pool_lm <- psfmi_lm(data=lbpmilr, nimp=5, impvar="Impnr", 
                      formula = Pain ~ Gender + Smoking + 
                      Function + JobControl + JobDemands + SocialSupport, 
                      keep.predictors = "Smoking", method="MPR", p.crit=0.05, 
                      direction="BW")
  
  pool_lm$RR_model_final
  pool_lm$multiparm_final
  pool_lm$predictors_out  

## -----------------------------------------------------------------------------

  library(psfmi)
  pool_lm <- psfmi_lm(data=lbpmilr, nimp=5, impvar="Impnr", 
                      formula = Pain ~ Gender + Smoking + 
                        Function + JobControl + factor(Carrying) + 
                        factor(Satisfaction) +
                        factor(Carrying):Smoking + Gender:Smoking, 
                      method="D2", p.crit=0.05, 
                      direction="BW")
  
  pool_lm$RR_model_final
  pool_lm$multiparm_final
  pool_lm$predictors_out 
  

## -----------------------------------------------------------------------------

  library(psfmi)
  pool_lm <- psfmi_lm(data=lbpmilr, nimp=5, impvar="Impnr", 
                      formula = Pain ~ Gender + Smoking + 
                      Function + JobControl + factor(Carrying) + factor(Satisfaction) +
                        factor(Carrying):Smoking + Gender:Smoking, 
                      keep.predictors = c("Smoking*Carrying", "JobControl"), method="D1", 
                      p.crit=0.05, direction="BW")
  
  pool_lm$RR_model_final
  pool_lm$multiparm_final
  pool_lm$predictors_out 
 

## -----------------------------------------------------------------------------

  library(psfmi)
  pool_lm <- psfmi_lm(data=lbpmilr, nimp=5, impvar="Impnr", 
                      formula = Pain ~ Gender + Smoking + 
                      JobControl + factor(Carrying) + factor(Satisfaction) +
                      factor(Carrying):Smoking + rcs(Function, 3), 
                      method="D1", 
                      p.crit=0.05, direction="BW")
  
  pool_lm$RR_model_final
  pool_lm$multiparm_final
  pool_lm$predictors_out 
  

