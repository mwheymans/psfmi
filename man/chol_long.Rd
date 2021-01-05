\name{chol_long}
\alias{chol_long}
\docType{data}
\title{Long dataset of persons from the The Amsterdam Growth and Health Longitudinal Study (AGHLS)}
\description{
Long dataset of persons from the The Amsterdam Growth and Health Longitudinal Study (AGHLS)
}
\usage{data(chol_long)}
\format{
  A data frame with 588 observations on the following 7 variables.
  \describe{
    \item{\code{ID}}{continuous}
    \item{\code{fitness}}{continuous}
    \item{\code{Smoking}}{dichotomous, 1=yes, 0=no}
    \item{\code{Sex}}{dichotomous}
    \item{\code{Time}}{categorical}
    \item{\code{Cholesterol}}{continuous}
    \item{\code{SumSkinfolds}}{continuous}
  }
}
\examples{
data(chol_long)
## maybe str(chol_long)
}
\keyword{datasets}
