---
title: "psfmi"
author: "Martijn W Heymans"
output: github_document
---

[![CRAN_Release_Badge](https://www.r-pkg.org/badges/version-ago/psfmi)](https://CRAN.R-project.org/package=psfmi)
[![Monthly downloads badge](https://cranlogs.r-pkg.org/badges/last-month/psfmi?color=blue)](https://CRAN.R-project.org/package=psfmi)
[![Travis Build Status](https://travis-ci.com/mwheymans/psfmi.svg?branch=master)](https://travis-ci.org/mwheymans/psfmi) 
[![minimal R version](https://img.shields.io/badge/R%3E%3D-3.6.0-6666ff.svg)](https://cran.r-project.org/)

The package provides functions to apply pooling, backward and forward selection 
of logistic, Cox regression and Multilevel (mixed models) prediction 
models in multiply imputed datasets. Backward and forward selection can be done 
from the pooled model using Rubin's Rules (RR), the D1, D2, D3 and 
the median p-values method. The model can contain	continuous, dichotomous, 
categorical predictors and interaction terms between all these type of predictors. 
Continuous predictors	can also be introduced as restricted cubic spline coefficients. 

The package includes a function to evaluate the stability 
of the models	using bootstrapping and cluster bootstrapping. The package further 
contains functions to generate pooled model performance measures in multiply 
imputed datasets as AUC value, R-squares, scaled Brier score, fit test values and 
calibration	plots for logistic regression models. 

Internal validation of the developed model can be performed with cross-validation or 
bootstrapping over multiply imputed datasets. The adjusted intercept after shrinkage 
of the pooled regression coefficients can be subsequently obtained. 
Backward and forward selection as part of internal validation is possible. 
Also a function to externally validate logistic	prediction models in 
multiple imputed datasets is available.

The package needs R 3.6.0 or higher. You can use the functions, 
after you have installed the package from CRAN or the Github website 
(development version). For Github you first have to install and activate 
the devtools package. Use the following code to install and activate 
the package from Github:

> install.packages("devtools")

> library(devtools)

> devtools::install_github("mwheymans/psfmi")

> library(psfmi)

Have fun! 





