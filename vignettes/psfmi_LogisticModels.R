## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------

  library(psfmi)
  pool_lr <- psfmi_lr(data=lbpmilr, nimp=5, impvar="Impnr", 
                      formula = Chronic ~ Gender + Smoking + 
                      Function + JobControl + JobDemands + SocialSupport, 
                      method="D1")
  
  pool_lr$RR_model
 

## -----------------------------------------------------------------------------

  library(psfmi)
  pool_lr <- psfmi_lr(data=lbpmilr, nimp=5, impvar="Impnr", 
                      formula = Chronic ~ Gender + Smoking + 
                      Function + JobControl + JobDemands + SocialSupport, 
                      keep.predictors = "Smoking", method="D1", p.crit=0.05, 
                      direction="BW")
  
  pool_lr$RR_model_final
  pool_lr$multiparm_final
  pool_lr$predictors_out
  

## -----------------------------------------------------------------------------

  library(psfmi)
  pool_lr <- psfmi_lr(data=lbpmilr, nimp=5, impvar="Impnr", 
                      formula = Chronic ~ Gender + Smoking + 
                      Function + JobControl + JobDemands + SocialSupport, 
                      keep.predictors = "Smoking", method="MPR", p.crit=0.05, 
                      direction="BW")
  
  pool_lr$RR_model_final
  pool_lr$multiparm_final
  pool_lr$predictors_out  
  

## -----------------------------------------------------------------------------

  library(psfmi)
  pool_lr <- psfmi_lr(data=lbpmilr, nimp=5, impvar="Impnr", 
                      formula = Chronic ~ Gender + Smoking + 
                        Function + JobControl + factor(Carrying) + 
                        factor(Satisfaction) +
                        factor(Carrying):Smoking + Gender:Smoking, 
                      method="D2", p.crit=0.05, 
                      direction="BW")
  
  pool_lr$RR_model_final
  pool_lr$multiparm_final
  pool_lr$predictors_out 
  

## -----------------------------------------------------------------------------

  library(psfmi)
  pool_lr <- psfmi_lr(data=lbpmilr, nimp=5, impvar="Impnr", 
                      formula = Chronic ~ Gender + Smoking + 
                      Function + JobControl + factor(Carrying) + factor(Satisfaction) +
                        factor(Carrying):Smoking + Gender:Smoking, 
                      keep.predictors = c("Smoking*Carrying", "JobControl"), method="D1", 
                      p.crit=0.05, direction="BW")
  
  pool_lr$RR_model_final
  pool_lr$multiparm_final
  pool_lr$predictors_out 
  

## -----------------------------------------------------------------------------

  library(psfmi)
  pool_lr <- psfmi_lr(data=lbpmilr, nimp=5, impvar="Impnr", 
                      formula = Chronic ~ Gender + Smoking + 
                      JobControl + factor(Carrying) + factor(Satisfaction) +
                      factor(Carrying):Smoking + rcs(Function, 3), 
                      method="D1", 
                      p.crit=0.05, direction="BW")
  
  pool_lr$RR_model_final
  pool_lr$multiparm_final
  pool_lr$predictors_out 


