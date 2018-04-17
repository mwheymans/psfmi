# psfmi

Welcome to the psfmi package. With this package you can pool 
logistic and Cox regression models after multiple imputation.
A function called psfmi_lr is available to pool logistic regression 
models and a function called psfmi_coxr is available to pool
right censored Cox regression models.

You can use these functions, after you have installed the
package from the the Github website. Before that you have to install 
and activate the devtools package. Use the following code for that:

install.packages("devtools")
library(devtools)
devtools::install_github("mwheymans/psfmi")
library(psfmi)
