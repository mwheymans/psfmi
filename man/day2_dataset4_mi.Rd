\name{day2_dataset4_mi}
\alias{day2_dataset4_mi}
\docType{data}
\title{Dataset of low back pain patients with missing values}
\description{
Dataset of low back pain patients with missing values in 2 variables
}
\usage{data(day2_dataset4_mi)}
\format{
  A data frame with 100 observations on the following 8 variables.
  \describe{
    \item{\code{ID}}{continuous: unique patient numbers}
    \item{\code{Pain}}{continuous: Pain intensity}
    \item{\code{Tampa}}{continuous: Fear of Movement scale}
    \item{\code{Function}}{continuous: Functional Status}
    \item{\code{JobSocial}}{continuous}
    \item{\code{FAB}}{continuous: Fear Avoidance Beliefs}
    \item{\code{Gender}}{dichotomous: 1 = male, 0 = female}
    \item{\code{Radiation}}{dichotomous: 1 = yes, 0 = no}
  }
}
\examples{
data(day2_dataset4_mi)
## maybe str(day2_dataset4_mi)
}
\keyword{datasets}
