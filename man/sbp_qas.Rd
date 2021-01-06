\name{sbp_qas}
\alias{sbp_qas}
\docType{data}
\title{Dataset with blood pressure measurements}
\description{
Dataset with blood pressure measurements
}
\usage{data(sbp_qas)}
\format{
  A data frame with 32 observations on the following 5 variables.
  \describe{
    \item{\code{pat_id}}{continuous}
    \item{\code{sbp}}{continuous: systolic blood pressure}
    \item{\code{bmi}}{continuous: body mass index}
    \item{\code{age}}{continuous: age (years)}
    \item{\code{smk}}{dichotomous: 0 = no, 1 = yes}
  }
}
\examples{
data(sbp_qas)
## maybe str(sbp_qas)
}
\keyword{datasets}
