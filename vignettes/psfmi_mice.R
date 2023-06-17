## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------

  library(psfmi)
  library(mice)

  imp <- mice(lbp_orig, m=5, maxit=5) 
  
  data_comp <- complete(imp, action = "long", include = FALSE)
  
  library(psfmi)
  pool_lr <- psfmi_lr(data=data_comp, nimp=5, impvar=".imp", 
                      formula=Chronic ~ Gender + Smoking + Function + 
                      JobControl + JobDemands + SocialSupport, method="D1")
  pool_lr$RR_model
  
 

## -----------------------------------------------------------------------------

  library(psfmi)
  library(mice)

  imp <- mice(lbp_orig, m=5, maxit=5) 
  
  data_comp <- complete(imp, action = "long", include = FALSE)
  
  library(psfmi)
  pool_lr <- psfmi_lr(data=data_comp, nimp=5, impvar=".imp", 
                      formula=Chronic ~ Gender + Smoking + Function + 
                      JobControl + JobDemands + SocialSupport, 
                      p.crit = 0.157, method="D1", direction = "FW")
  
  pool_lr$RR_model_final
 

