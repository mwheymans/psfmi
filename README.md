
<!-- README.md is generated from README.Rmd. Please edit that file -->

# psfmi

[![CRAN_Release_Badge](https://www.r-pkg.org/badges/version-ago/psfmi)](https://CRAN.R-project.org/package=psfmi)
[![Monthly downloads
badge](https://cranlogs.r-pkg.org/badges/last-month/psfmi?color=blue)](https://CRAN.R-project.org/package=psfmi)
[![Travis Build
Status](https://travis-ci.com/mwheymans/psfmi.svg?branch=master)](https://travis-ci.org/mwheymans/psfmi)
[![minimal R
version](https://img.shields.io/badge/R%3E%3D-4.0.0-6666ff.svg)](https://cran.r-project.org/)

The package provides functions to apply pooling, backward and forward
selection of linear, logistic and Cox regression models across multiply
imputed data sets using Rubin’s Rules (RR). The D1, D2, D3, D4 and the
median p-values method can be used to pool the significance of
categorical variables (multiparameter test). The model can contain
continuous, dichotomous, categorical and restricted cubic spline
predictors and interaction terms between all these type of variables.
Variables can also be forced in the model during selection.

Validation of the prediction models can be performed with
cross-validation or bootstrapping across multiply imputed data sets and
pooled model performance measures as AUC value, R-square, Hosmer and
Lemeshow test, scaled Brier score and calibration plots are generated.
Also a function to externally validate logistic prediction models across
multiple imputed data sets is available and a function to compare models
in multiply imputed data.

## Installation

You can install the released version of psfmi with:

``` r
install.packages("psfmi")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mwheymans/psfmi")
```

## Citation

Cite the package as:

``` r

Martijn W Heymans (2021). psfmi: Prediction Model Pooling, Selection and Performance Evaluation 
Across Multiply Imputed Datasets. R package version 1.1.0. https://mwheymans.github.io/psfmi/ 
```

## Examples

This example shows you how to pool a logistic regression model across 5
multiply imputed datasets and that includes two restricted cubic spline
variables and a categorical, continuous and dichotomous variable. The
pooling method that is used is method D1.

``` r
library(psfmi)

pool_lr <- psfmi_lr(data=lbpmilr, formula = Chronic ~ rcs(Pain, 3) + 
                      JobDemands + rcs(Tampascale, 3) + factor(Satisfaction) + 
                      Smoking, nimp=5, impvar="Impnr", method="D1")

pool_lr$RR_model
#> $`Step 1 - no variables removed -`
#>                            term      estimate  std.error  statistic        df
#> 1                   (Intercept) -21.374498123 7.96491209 -2.6835824  65.71094
#> 2                    JobDemands  -0.007500147 0.05525835 -0.1357288  38.94021
#> 3                       Smoking   0.072207184 0.51097303  0.1413131  47.98415
#> 4         factor(Satisfaction)2  -0.506544055 0.56499941 -0.8965391 139.35335
#> 5         factor(Satisfaction)3  -2.580503376 0.77963853 -3.3098715 100.66273
#> 6              rcs(Pain, 3)Pain  -0.090675006 0.50510774 -0.1795162  26.92182
#> 7             rcs(Pain, 3)Pain'   1.183787048 0.55697046  2.1254036  94.79276
#> 8  rcs(Tampascale, 3)Tampascale   0.583697990 0.22707747  2.5704796  77.83368
#> 9 rcs(Tampascale, 3)Tampascale'  -0.602128298 0.29484065 -2.0422160  31.45559
#>       p.value           OR    lower.EXP  upper.EXP
#> 1 0.009206677 5.214029e-10 6.460344e-17 0.00420815
#> 2 0.892734942 9.925279e-01 8.875626e-01 1.10990663
#> 3 0.888214212 1.074878e+00 3.847422e-01 3.00295282
#> 4 0.371511077 6.025744e-01 1.971829e-01 1.84141687
#> 5 0.001296125 7.573587e-02 1.612863e-02 0.35563604
#> 6 0.858876729 9.133145e-01 3.239353e-01 2.57503035
#> 7 0.036152843 3.266722e+00 1.081155e+00 9.87043962
#> 8 0.012063538 1.792655e+00 1.140659e+00 2.81733025
#> 9 0.049589266 5.476448e-01 3.002599e-01 0.99885104

pool_lr$multiparm
#> $`Step 1 - no variables removed -`
#>                      p-values D1 F-statistic
#> JobDemands           0.892487763  0.01842230
#> Smoking              0.887968553  0.01996939
#> factor(Satisfaction) 0.002611518  6.04422205
#> rcs(Pain,3)          0.014630986  4.84409246
#> rcs(Tampascale,3)    0.130741167  2.24870192
```

This example shows you how to apply forward selection of the above model
using a p-value of 0.05.

``` r
library(psfmi)

pool_lr <- psfmi_lr(data=lbpmilr, formula = Chronic ~ rcs(Pain, 3) + 
                      JobDemands + rcs(Tampascale, 3) + factor(Satisfaction) + 
                      Smoking, p.crit = 0.05, direction="FW", 
                      nimp=5, impvar="Impnr", method="D1")
#> Entered at Step 1 is - rcs(Pain,3)
#> Entered at Step 2 is - factor(Satisfaction)
#> 
#> Selection correctly terminated, 
#> No new variables entered the model

pool_lr$RR_model_final
#> $`Final model`
#>                    term   estimate std.error  statistic        df     p.value
#> 1           (Intercept) -3.6027668 1.5427414 -2.3353018  60.25659 0.022875170
#> 2 factor(Satisfaction)2 -0.4725289 0.5164342 -0.9149838 145.03888 0.361718841
#> 3 factor(Satisfaction)3 -2.3328994 0.7317131 -3.1882707 122.95905 0.001815476
#> 4      rcs(Pain, 3)Pain  0.6514983 0.4028728  1.6171315  51.09308 0.112008088
#> 5     rcs(Pain, 3)Pain'  0.4703811 0.4596490  1.0233483  75.29317 0.309419924
#>           OR   lower.EXP upper.EXP
#> 1 0.02724823 0.001245225 0.5962503
#> 2 0.62342367 0.224644070 1.7301016
#> 3 0.09701406 0.022793375 0.4129150
#> 4 1.91841309 0.854476033 4.3070942
#> 5 1.60060402 0.640677978 3.9987846

pool_lr$multiparm
#> $`Step 0 - selected - rcs(Pain,3)`
#>                        p-value D1
#> JobDemands           7.777737e-01
#> Smoking              9.371529e-01
#> factor(Satisfaction) 9.271071e-01
#> rcs(Pain,3)          3.282999e-07
#> rcs(Tampascale,3)    2.780012e-06
#> 
#> $`Step 1 - selected - factor(Satisfaction)`
#>                       p-value D1
#> JobDemands           0.952900908
#> Smoking              0.769394518
#> factor(Satisfaction) 0.004738608
#> rcs(Tampascale,3)    0.125280292
```

More examples for logistic, linear and Cox regression models as well as
internal and external validation of prediction models can be found on
the [package website](https://mwheymans.github.io/psfmi/) or in the
online book [Applied Missing Data
Analysis](https://bookdown.org/mwheymans/bookmi/).
