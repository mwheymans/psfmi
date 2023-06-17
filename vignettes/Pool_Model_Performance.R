## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------

library(psfmi)

perf <- pool_performance(data=lbpmilr, nimp=5, impvar="Impnr", 
  formula = Chronic ~ Gender + Pain + Tampascale + 
  Smoking + Function + Radiation + Age + 
  Duration + BMI, 
  cal.plot=TRUE, plot.method="mean", 
  groups_cal=10, model_type="binomial")
  
perf


