# Package 1.4.0

* Added a new function nri_cox.R for reclassification analyses of Cox models.
* Added the possibility to pool stratified Cox models in the function psfmi_coxr.

# Package 1.1.0

* Added df_com to the D1 and D2 method for small sample
correction.
* Update strata option in functions cv_MI and cv_MI_RR.

# Package 1.0.0

* Added the functions psfmi_lm, psfmi_lm_fw, psfmi_lr_bw for pooling
and backward and forward selection of linear regression models.
* Added the functions glm_bw, glm_fw, coxph_bw and coxph_fw for 
backward and forward selection of linear, logistic and Cox models 
in a single dataset based on the likelihood ratio statistic.
* Function psfmi_perform is now deprecated, use psfmi_validate instead.
* Function bw_single is now deprecated, use glm_bw instead.
* Added the function hoslem_test and implemented this in the 
function pool_performance.
* added pooled concordance and R-squared measures for Cox regression to
function pool_performance.
* Added the function pool_D2, to pool chi-square statistics.
* Added the function pool_D4, to pool likelihood ratio tests.
* Added the internal function pool_performance_internal, used internally 
by psfmi_perform.
* Option plot.indiv in function pool_performance and mivalext_lr is deprecated, 
use plot.method instead.
* created a new vignette for the psfmi_lm function and updated
the other ones.
* corrected other bug fixes.
* Updated package website with pkgdown (included citation).

# Package 0.7.1

* Added the functions pool_compare_models and pool_reclassification.
* applied some bug fixes in the function bw_single where selection with anova
was replaced with Anova.
* added internal function RR_diff_prop to pool difference in proportions and
related SE with RR.
* Updated vignettes.
* Updated  package website with pkgdown.
* Added datasets to be used in tutorials.

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
