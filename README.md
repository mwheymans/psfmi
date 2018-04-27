# psfmi

Welcome to the psfmi package. 

With this package you can pool logistic or Cox regression models,
or perform backward variable selection, in multiply imputed datasets. 
The models may include continuous, dichotomous and categorical 
predictors. Also interaction terms between different type of 
predictor variables are allowed. It is also possible to force  
predictors or interaction terms in the model during predictor selection.

The pooling methods are the Meng and Rubin likelihood ratio statistics 
method, the pooling covariance matrix or D1 method, the pooling of 
Chi-square values or D2 method and a method that is called 
Median P Rule (MPR) that pools the median of the p-values.

A function called psfmi_lr is available for logistic regression 
models and another function that is called psfmi_coxr, 
for right censored Cox regression models.

A function that is called miperform_lr evaluates the apparent performance 
of logistic regression prediction models in imputed datasets. 
The performanance measures that are reported are the ROC/AUC, 
(Nagerkerke) R-squared and the Hosmer and Lemeshow test. 
Also calibration curves can be generated, overlayed curves that pool 
all calibration curves that are estimated in each imputed datset, 
or individual curves as a result of calibrating the model in each 
imputed dataset separately. The pooled linear predictor is also returned.

You can use these functions, after you have installed the
package from the Github website. Before that you have to install 
and activate the devtools package. Use the following code to
install and activate the package:

> install.packages("devtools")

> library(devtools)

> devtools::install_github("mwheymans/psfmi")

> library(psfmi)

