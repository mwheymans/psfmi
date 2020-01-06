\name{ipdna_md}
\alias{ipdna_md}
\docType{data}
\title{Example dataset for the psfmi_mm function}
\description{
5 imputed datasets of the first 10 centres of the IPDNa dataset in the micemd package.
}
\usage{data(ipdna_md)}
\format{
  A data frame with 13390 observations on the following 13 variables.
  \describe{
    \item{\code{.imp}}{a numeric vector}
    \item{\code{.id}}{a numeric vector}
    \item{\code{centre}}{cluster variable}
    \item{\code{gender}}{dichotomous}
    \item{\code{bmi}}{continuous}
    \item{\code{age}}{continuous}
    \item{\code{sbp}}{continuous}
    \item{\code{dbp}}{continuous}
    \item{\code{hr}}{continuous}
    \item{\code{lvef}}{dichotomous}
    \item{\code{bnp}}{categorical}
    \item{\code{afib}}{continuous}
    \item{\code{bmi_cat}}{categorical}
  }
}
\examples{
data(ipdna_md)
## maybe str(ipdna_md)

#summary per study
by(ipdna_md, ipdna_md$centre, summary)
}
\keyword{datasets}
