---
title: "cran_comments"
---
## Resubmission
This is a resubmission. In this version I have:
* Updated version nr to 1.3.0.
* Added a new function nri_cox.R for reclassification analyses of Cox models.
* Added the possibility to pool stratified Cox models in the function psfmi_coxr.

── R CMD check results ────────────────────────────────────────────── psfmi 1.3.0 ────
Duration: 3m 16.3s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## Test environments
*	R studio 2023.06.0, R version 4.3.0
* win-builder Development: * DONE Status: OK
* win-builder Release: * DONE Status: OK

## Resubmission
This is a resubmission. In this version I have:
* Updated version nr to 1.1.0.
* Added df_com to the D1 and D2 method for small sample
correction.
* Defined strata=strata in functions cv_MI and cv_MI_RR.

## devtools::test()
ℹ Testing psfmi
✔ | F W S  OK | Context
⠏ |         0 | testthat.R                                                                                    
══ Results ═══════════════════════════════════════════════════════════════════════════════════════════════════
Duration: 0.2 s

[ FAIL 0 | WARN 0 | SKIP 0 | PASS 0 ]

Nice code.

## R CMD check results
── R CMD check results ──────────────────────────────────────────────────────────────────────── psfmi 1.1.0 ────
Duration: 5m 22.6s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## Test environments
*	R studio 2022.07.1, R version 4.2.1
* win-builder Development: * DONE Status: OK
* win-builder Release: * DONE Status: OK

## Resubmission
This is a resubmission. In this version I have:
* Updated version nr to 1.0.0.
* Added the functions psfmi_lm, psfmi_lm_fw, psfmi_lr_bw for pooling
and backward and forward selection of linear regression models.
* Added the functions glm_bw, glm_fw, coxph_bw and coxph_fw for 
backward and forward selection of linear, logistic and Cox models 
in a single dataset based on the likelihood ratio statistic.
* Function psfmi_perform is now deprecated, use psfmi_validate instead.
* Function bw_single is now deprecated, use glm_bw instead.
* Added the function hoslem_test and implemented this in the 
function pool_performance.
* Added pooled concordance and R-squared measures for Cox regression to
function pool_performance.
* Added the function pool_D2, to pool chi-square statistics.
* Added the function pool_D4, to pool likelihood ratio tests.
* Added the internal function pool_performance_internal, used internally 
by psfmi_perform.
* Option plot.indiv in function pool_performance and mivalext_lr is deprecated, 
use plot.method instead.
* Created a new vignette for the psfmi_lm function and updated
the other ones.
* corrected other bug fixes.
* Updated package website with pkgdown (included citation).

## R CMD check results
Duration: 6m 8.1s

0 errors √ | 0 warnings √ | 0 notes √

## devtools::test()
i Loading psfmi
i Testing psfmi
√ |  OK F W S | Context
/ |   0       | testthat.R                                                                                  
== Results =================================================================================================
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 0 ]

Nice code.

## Test environments
*	R studio 1.4.1717, R version 4.1.1
* win-builder Development: * DONE Status: OK
* win-builder Release: * DONE Status: OK
* on Travis CI: R session information
Running under: Ubuntu 16.04.6 LTS
Done. Your build exited with 0.
  

## Resubmission
This is a resubmission. In this version I have:
Updated version nr to 0.7.1.
Added new functions to compare models, the functions pool_compare_models and 
pool_reclassification. Applied some bug fixes in the function bw_single. 
Added internal function RR_diff_prop to pool difference in proportions. 
Updated vignettes and package website with pkgdown and added datasets to 
be used in tutorials.

## R CMD check results
Duration: 5m 19.8s

0 errors √ | 0 warnings √ | 0 notes √

## devtools::test()
Loading psfmi
Testing psfmi
√ |  OK F W S | Context
/ |   0       | testthat.R                                                                                              
== Results =============================================================================================================
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 0 ]

## Test environments
*	R studio 1.4.970, R version 4.0.3
* win-builder Development: * DONE Status: OK
* win-builder Release: * DONE Status: OK
* on Travis CI: R session information
$ Rscript -e 'sessionInfo()'
R version 4.0.2 (2020-06-22)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 16.04.6 LTS
Done. Your build exited with 0.
  
## Resubmission
This is a resubmission. In this version I have:
Updated version nr to 0.5.0.
Added new methods to combine cross-validation with multiple imputation
in function psfmi_perform and added the possibility for forward 
and backward selection in the functions psfmi_lr and psfmi_coxr. 
Updated the vignettes. 

## R CMD check results
Duration: 4m 55s

0 errors √ | 0 warnings √ | 1 note x

checking installed package size ... NOTE
    installed size is  7.0Mb
    sub-directories of 1Mb or more:
      doc   6.5Mb

Known note.

## devtools::test()
Loading psfmi
Testing psfmi
√ |  OK F W S | Context
√ |   0       | testthat.R

== Results =====================================================================
OK:       0
Failed:   0
Warnings: 0
Skipped:  0

Nice code.

## Test environments
*	R studio 1.3.1073, R version 4.0.2
* win-builder Development: Status: OK, R version 4.0.2
* win-builder Release: Status: OK, R version 4.0.2
* on Travis CI: Operating System Details Distributor ID:	Ubuntu
  Description:	Ubuntu 16.04.6 LTS
  Release:	16.04
  Codename:	xenial (on Travis CI) (Done. Your build exited with 0)

## Resubmission
This is a resubmission. In this version I have:
Updated version nr to 0.2.0.
Added the new functions pool_intadj.R, pool_performance.R, psfmi_mm.R,
psfmi_mm_multiparm, psfmi_perform, psfmi_stab.
Excluded function miperform_lr. 
Updated the vignettes. 

## R CMD check results
Duration: 2m 45.6s

0 errors √ | 0 warnings √ | 0 notes √

## devtools::test()
Loading psfmi
Testing psfmi
√ |  OK F W S | Context
√ |   0       | testthat.R

== Results =====================================================================
OK:       0
Failed:   0
Warnings: 0
Skipped:  0

## Resubmission
This is a resubmission. In this version I have:
Added small executable examples in all the exported functions'
Rd files to illustrate the use of the exported function.
Deleted comments in examples.

## Resubmission
This is a resubmission. In this version I have:
Provided small executable examples. 
Omitted cat()/print() in functions and used message() on a minority basis.
Changed the functions so that the user can extract separate objects. 
## Test environments
*	MacOS Mojave versie 10.14.4, R studio 1.1.383, R version 3.5.3
* win-builder (devel and release), R version 3.6.0
* Ubuntu 14.04.5 LTS (on Travis CI) (Done. Your build exited with 0)

## R CMD check results
Duration: 1m 29.8s

0 errors √ | 0 warnings √ | 0 notes √

## devtools::test()
Loading psfmi
Testing psfmi
√ | OK F W S | Context
√ |  0       | testthat.R
== Results =====================================================================
OK:       0
Failed:   0
Warnings: 0
Skipped:  0

## Resubmission
This is a resubmission. In this version I have:
Excluded annotated lines at location of example.

## Resubmission
This is a resubmission. In this version I have:
Started the Description with just "Provides..." and did not capitalize non-names (e.g. Likelihood - likelihood).

I have added references as  (year) <doi:...>  

I have changed my examples into small executable examples in all of my Rd-files.

## Resubmission
This is a resubmission. In this version I have:
Omitted the mis-spelled words, the date field and adjusted the examples so that they run <5s.

## Test environments
*	MacOS Mojave versie 10.14.4, R studio 1.1.383, R version 3.5.3
* win-builder (devel and release)
* Ubuntu 14.04.5 LTS (on Travis CI) (Done. Your build exited with 0)

## R CMD check results
errors √ | 0 warnings √ | 0 notes √

## devtools::test()
Loading psfmi
Testing psfmi
√ | OK F W S | Context
√ |  0       | testthat.R
== Results =====================================================================
OK:       0
Failed:   0
Warnings: 0
Skipped:  0

You are a coding rockstar!
