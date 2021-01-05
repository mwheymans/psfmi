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
    \item{\code{sbp}}{continuous: systo¬lic blood pressure (SBP, in mm Hg)}
    \item{\code{bmi}}{continuous: body mass index (BMI = weight / height2, weight in kg, height in m)}
    \item{\code{age}}{continuous: age (AGE, in years)}
    \item{\code{smk}}{dichotomous: SMK; 0 = no, 1 = yes}
  }
}
\examples{
data(sbp_qas)
## maybe str(sbp_qas)
}
\keyword{datasets}
