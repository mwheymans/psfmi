% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coxph_fw.R
\name{coxph_fw}
\alias{coxph_fw}
\title{Predictor selection function for forward selection of Cox regression models in 
 single complete dataset.}
\usage{
coxph_fw(
  data,
  formula = NULL,
  status = NULL,
  time = NULL,
  predictors = NULL,
  p.crit = 1,
  cat.predictors = NULL,
  spline.predictors = NULL,
  int.predictors = NULL,
  keep.predictors = NULL,
  nknots = NULL
)
}
\arguments{
\item{data}{A data frame.}

\item{formula}{A formula object to specify the model as normally used by coxph.
See under "Details" and "Examples" how these can be specified.}

\item{status}{The status variable, normally 0=censoring, 1=event.}

\item{time}{Survival time.}

\item{predictors}{Character vector with the names of the predictor variables.
At least one predictor variable has to be defined. Give predictors unique names
and do not use predictor name combinations with numbers as, age2, gnder10, etc.}

\item{p.crit}{A numerical scalar. P-value selection criterium. A value of 1
provides the pooled model without selection.}

\item{cat.predictors}{A single string or a vector of strings to define the
categorical variables. Default is NULL categorical predictors.}

\item{spline.predictors}{A single string or a vector of strings to define the
(restricted cubic) spline variables. Default is NULL spline predictors. See details.}

\item{int.predictors}{A single string or a vector of strings with the names of the variables that form
an interaction pair, separated by a “:” symbol.}

\item{keep.predictors}{A single string or a vector of strings including the variables that are forced
in the model during predictor selection. All type of variables are allowed.}

\item{nknots}{A numerical vector that defines the number of knots for each spline predictor separately.}
}
\value{
An object of class \code{smods} (single models) from
 which the following objects can be extracted: original dataset as \code{data}, final selected
 model as \code{RR_model_final}, model at each selection step \code{RR_model},
 p-values at final step \code{multiparm_final}, and at each step as \code{multiparm}, 
 formula object at final step as \code{formula_final}, 
 and at each step as \code{formula_step} and for start model as \code{formula_initial}, 
 predictors included at each selection step as \code{predictors_in}, predictors excluded
 at each step as \code{predictors_out}, and \code{time}, \code{status}, \code{p.crit}, \code{call},
 \code{model_type}, \code{predictors_final} for names of predictors in final selection step and 
 \code{predictors_initial} for names of predictors in start model and \code{keep.predictors} for
 variables that are forced in the model during selection.
}
\description{
\code{coxph_bw} Forward selection of Cox regression models in single complete 
 dataset using as selection method the partial likelihood-ratio statistic.
}
\details{
A typical formula object has the form \code{Surv(time, status) ~ terms}. Categorical variables has to
 be defined as \code{Surv(time, status) ~ factor(variable)}, restricted cubic spline variables as
 \code{Surv(time, status) ~ rcs(variable, 3)}. Interaction terms can be defined as
 \code{Surv(time, status) ~ variable1*variable2} or \code{Surv(time, status) ~ variable1 + variable2 + 
 variable1:variable2}. All variables in the terms part have to be separated by a "+".
}
\examples{
lbpmicox1 <- subset(psfmi::lbpmicox, Impnr==1) # extract first imputed dataset
res_single <- coxph_bw(data=lbpmicox1, p.crit = 0.05, formula=Surv(Time, Status) ~
                           Previous +  Radiation + Onset + Age + Tampascale + 
                           Pain + JobControl + factor(Satisfaction), 
                           spline.predictors = "Function",
                           nknots = 3)
         
res_single$RR_model_final
res_single$multiparm_final

}
\references{
http://missingdatasolutions.rbind.io/
}
\author{
Martijn Heymans, 2021
}
