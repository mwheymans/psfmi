% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stab_single.R
\name{stab_single}
\alias{stab_single}
\title{Function to evaluate bootstrap predictor and model stability.}
\usage{
stab_single(pobj, nboot = 20, p.crit = 0.05, start_model = TRUE)
}
\arguments{
\item{pobj}{An object of class \code{smods} (single models), produced by a previous call to
\code{glm_bw}.}

\item{nboot}{A numerical scalar. Number of bootstrap samples to evaluate the stability. Default is 20.}

\item{p.crit}{A numerical scalar. Used as P-value selection criterium during bootstrap model selection.}

\item{start_model}{If TRUE the bootstrap evaluation takes place from the start model of object pobj, if
FALSE the final model is used for the evaluation.}
}
\value{
A \code{psfmi_stab} object from which the following objects can be extracted: bootstrap
 inclusion (selection) frequency of each predictor \code{bif}, total number each predictor is
 included in the bootstrap samples as \code{bif_total}, percentage a predictor is selected
 in each bootstrap sample as \code{bif_perc} and number of times a prediction model is selected in
 the bootstrap samples as \code{model_stab}.
}
\description{
\code{stab_single} Stability analysis of predictors and prediction models selected with
 the \code{glm_bw}.
}
\details{
The function evaluates predictor selection frequency in bootstrap samples.
 It uses as input an object of class \code{smods} as a result of a
 previous call to the \code{glm_bw}.
}
\examples{
 model_lr <- glm_bw(formula = Radiation ~ Pain + factor(Satisfaction) + 
   rcs(Tampascale,3) + Age + Duration + JobControl + JobDemands + SocialSupport, 
   data=lbpmilr_dev, p.crit = 0.05)

\dontrun{
 stab_res <- stab_single(model_lr, start_model = TRUE, nboot=20, p.crit=0.05)
 stab_res$bif
 stab_res$bif_perc
 stab_res$model_stab
}

}
\references{
Heymans MW, van Buuren S. et al. Variable selection under multiple imputation using the bootstrap
  in a prognostic study. BMC Med Res Methodol. 2007;13:7-33.

Sauerbrei W, Schumacher M. A bootstrap resampling procedure for model building:
  application to the Cox regression model. Stat Med. 1992;11:2093–109.

Royston P, Sauerbrei W (2008) Multivariable model-building – a pragmatic approach to
  regression analysis based on fractional polynomials for modelling continuous variables. (2008).
  Chapter 8, Model Stability. Wiley, Chichester.

Heinze G, Wallisch C, Dunkler D. Variable selection - A review and
 recommendations for the practicing statistician. Biom J. 2018;60(3):431-449.

http://missingdatasolutions.rbind.io/
}
