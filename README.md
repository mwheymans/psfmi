# psfmi

The psfmi package 

With this package you can pool logistic or Cox regression models and
(generalized) linear mixed models in multiply imputed datasets
or perform backward variable selection from the pooled model. 
The models may include continuous, dichotomous, categorical (> 2 
categories) and spline predictors. Also interaction terms between these type of 
predictor variables are possible. It is also possible to force (spline)  
predictors or interaction terms in the model during predictor selection.
This is the first package that contains these possibilities. 

A function called psfmi_lr is available for logistic regression 
models, a function that is called psfmi_coxr, for right censored 
Cox regression models and a function that is called psfmi_mm for
linear and logistic multilevel models (mixed models).

The basic pooling method is Rubin's Rules (RR). For categorical and 
spline predictors the following pooling methods are available: 
D1, D2, and D3 and a method that is called Median P Rule (MPR) 
that pools the median of p-values. The latter is a simple 
but very promising procedure.

With respect to introducing interaction terms in the model, only 
two-way interactions are allowed. If interaction terms are included 
and backward selection is applied, interaction terms are dropped 
from the model following the hierarchy principle. This means
that interactions are considered during backward selection when both
main effects are in the model.

The package also contains a function to study the stability of the
models during backward selection. This function is called psfmi_stab
and uses single bootstrapping for logistic and Cox regression
models and cluster bootstrapping for multilevel models.

A function that is called psfmi_perform evaluates the pooled performance 
of logistic regression prediction models in multiple imputed datasets. 
The performanance measures that are reported are the pooled ROC/AUC, 
(Nagerkerke) R-squared and the Brier score. Also calibration curves 
can be generated, overlayed curves that pool all calibration curves that 
are estimated in each imputed dataset, or individual curves as a result 
of calibrating the model in each imputed dataset separately. 
The pooled linear predictor is also returned. With this function it is also 
possible to apply internal validation using bootstrapping. Two methods
are available: One method, boot_MI, first draws	bootstrap samples and	
subsequently performs multiple imputation and with the other method, MI_boot, 
first bootstrap samples are drawn from each imputed dataset before results 
are combined. Backward selection as part of internal validation is also an option.

A function with the name mivalext_lr can be used to externally validate
prediction models in multiple imputed datasets. The following information 
of the externally validated model is provided: pooled ROC/AUC, (Nagelkerke) 
R-Square value, Hosmer and Lemeshow Test, pooled coefficients when the model 
is freely estimated in imputed datasets and the pooled linear predictor (LP), 
with information about miscalibration in intercept and slope. 

The package needs R 3.5.0 or higher. You can use the functions, 
after you have installed the package from CRAN or the Github website 
(development version). For Github you first have to install and activate 
the devtools package. Use the following code to install and activate 
the package from Github:

> install.packages("devtools")

> library(devtools)

> devtools::install_github("mwheymans/psfmi")

> library(psfmi)

