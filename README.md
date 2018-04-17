# psfmi

Welcome to the psfmi package. 

With this package you can pool logistic or Cox regression models,
or do backward variable selection with these models, 
in multiply imputed datasets. These models may include 
continuous, dichotomous and categorical predictors. Also
interaction terms between these type of variables are allowed.
It is also possible to force specific predictors or interaction 
terms in the model during predictor selection.

The pooling methods are the Meng and Rubin likelihood ratio statistics 
method, pooling covariance matrix or D1 method, pooling of 
Chi-square values or D2 method and a method that is called Median P Rule
that pools by taking the median of the p-values.

A function called psfmi_lr is available for logistic regression 
models and another function that is called psfmi_coxr, 
for right censored Cox regression models.

You can use these functions, after you have installed the
package from the Github website. Before that you have to install 
and activate the devtools package. Use the following code to
install and activate the package:

> install.packages("devtools")
> library(devtools)
> devtools::install_github("mwheymans/psfmi")
> library(psfmi)
