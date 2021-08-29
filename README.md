
<!-- README.md is generated from README.Rmd. Please edit that file -->

# psfmi

[![CRAN\_Release\_Badge](https://www.r-pkg.org/badges/version-ago/psfmi)](https://CRAN.R-project.org/package=psfmi)
[![Monthly downloads
badge](https://cranlogs.r-pkg.org/badges/last-month/psfmi?color=blue)](https://CRAN.R-project.org/package=psfmi)
[![Travis Build
Status](https://travis-ci.com/mwheymans/psfmi.svg?branch=master)](https://travis-ci.org/mwheymans/psfmi)
[![minimal R
version](https://img.shields.io/badge/R%3E%3D-4.0.0-6666ff.svg)](https://cran.r-project.org/)

The package provides functions to apply pooling, backward and forward
selection of linear, logistic and Cox regression prediction models in
multiply imputed data sets using Rubin’s Rules (RR), the D1, D2, D3 and
the median p-values method. The model can contain continuous,
dichotomous, categorical and restricted cubic spline predictors and
interaction terms between all these type of predictors.

Validation of the prediction models can be performed with
cross-validation or bootstrapping in multiply imputed data sets and
pooled model performance measures as AUC value, R-square, scaled Brier
score and calibration plots are generated. Also a function to externally
validate logistic prediction models in multiple imputed data sets is
available and a function to compare models in multiply imputed data.

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
Martijn Heymans (2021). psfmi: Prediction Model Selection and Performance Evaluation in
Multiple Imputed Datasets. R package version 1.0.0. https://mwheymans.github.io/psfmi/
```

## Example

This example shows you how to apply forward selection with a model that
includes a restricted cubic spline function and an interaction term.

``` r
library(psfmi)

pool_lr <- psfmi_lr(data=lbpmilr, formula = Chronic ~ rcs(Pain, 3) + JobDemands + rcs(Tampascale, 3) +
                   factor(Satisfaction) + Smoking + factor(Satisfaction)*rcs(Pain, 3) ,
                   p.crit = 0.05, direction="FW", nimp=5, impvar="Impnr",
                   method="D1")
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
```
