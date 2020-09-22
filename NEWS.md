# Package 0.5.0

* Added forward selection to the psfmi_lr and psfmi_coxr functions.
* Added the functions psfmi_lr_fw, psfmi_lr_bw, psfmi_coxr_fw and 
  psfmi_coxr_bw that are called by psfmi_lr or psfmi_coxr for 
  forward and backward selection after MI with logistic and Cox 
  regression models.
* Extended the function psfmi_perform by including the possibility 
  of using cross-validation in combination with multiple imputation.
* Added the functions cv_MI, cv_MI_RR, MI_cv_naive, boot_MI and MI_boot,
  that are called by psfmi_perform to combine cross-validation 
  or bootstrapping with MI.  
* Added the function bw_single for backward selection in a single
  dataset.
* Added the functions pool_auc, rsq_nagel and scaled_brier.
* Added the internal functions mean_auc_log, clean_P, miceImp.
* Updated vignettes.
* Created a package website with pkgdown.

# Package 0.2.0

* Added functions pool_intadj.R, pool_performance.R, psfmi_mm.R,
 psfmi_mm_multiparm, psfmi_perform, psfmi_stab.
* Added the class statement smodsmi to the functions psfmi_lr, psfmi_coxr  
 and psfmi_mm.
* Added the following output object information to psfmi_lr, psfmi_coxr  
 and psfmi_mm: model_type, predictors_in, predictors_out, fit.formula, 
 predictors_final and predictors_initial.
* Excluded function miperform_lr.
* Changed name object multiparm_p into multiparm.
* Updated vignettes. 
