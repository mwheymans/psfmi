\name{lbpmicox}
\alias{lbpmicox}
\docType{data}
\title{Example dataset for psfmi_coxr function}
\description{
10 imputed datasets
}
\usage{data(lbpmicox)}
\format{
  A data frame with 2650 observations on the following 18 variables.
  \describe{
    \item{\code{Impnr}}{a numeric vector}
    \item{\code{patnr}}{a numeric vector}
    \item{\code{Status}}{dichotomous event }
    \item{\code{Time}}{continuous follow up time variable}
    \item{\code{Duration}}{continuous}
    \item{\code{Previous}}{dichotomous}
    \item{\code{Radiation}}{dichotomous}
    \item{\code{Onset}}{dichotomous}
    \item{\code{Age}}{continuous}
    \item{\code{Tampascale}}{continuous}
    \item{\code{Pain}}{continuous}
    \item{\code{Function}}{continuous}
    \item{\code{Satisfaction}}{categorical}
    \item{\code{JobControl}}{continuous}
    \item{\code{JobDemand}}{continuous}
    \item{\code{Social}}{continuous}
    \item{\code{Expectation}}{a numeric vector}
    \item{\code{Expect_cat}}{categorical}
  }
}
\examples{
data(lbpmicox)
## maybe str(lbpmicox)
}
\keyword{datasets}
