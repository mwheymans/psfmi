\name{aortadis}
\alias{aortadis}
\docType{data}
\title{Dataset of patients with a aortadissection}
\description{
Original dataset of patients with a aortadissection
}
\usage{data(aortadis)}
\format{
  A data frame with 226 observations on the following 10 variables.
  \describe{
    \item{\code{Gender}}{dichotomous, 1=yes, 0=no}
    \item{\code{Age}}{continuous}
    \item{\code{Age_C}}{categorical: 0 = < 50 years, 1 = 50-59 years,
		2 = 60-69 years, 3 = 70-79 years, 4 = 80 years and older}
    \item{\code{Aortadis}}{dichotomous, 1=yes, 0=no}
    \item{\code{Acute}}{dichotomous, 1=yes, 0=no}
    \item{\code{Acute3}}{categorical: 0 = No, 1 = Little, 2 = Much}
    \item{\code{Stomach_Ache}}{dichotomous, 1=yes, 0=no}
    \item{\code{Smoking}}{dichotomous, 1=yes, 0=no}
    \item{\code{Radiation}}{dichotomous, 1=yes, 0=no}
  }
}
\examples{
data(aortadis)
## maybe str(aortadis)
}
\keyword{datasets}
