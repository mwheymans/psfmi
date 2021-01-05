\name{anderson}
\alias{anderson}
\docType{data}
\title{Data from a placebo-controlled RCT with leukemia patients}
\description{
Data from a placebo-controlled RCT with leukemia patients
}
\usage{data(anderson)}
\format{
  A data frame with 348 observations on the following 5 variables.
  \describe{
    \item{\code{remission}}{continuous:remission in weeks}
    \item{\code{status}}{dichotomous}
    \item{\code{treatment}}{dichotomous: 0=placebo, 1=verum}
    \item{\code{sex}}{dichotomous: 0=female, 1=male}
    \item{\code{log_wbc}}{continuous: Log (number of white blood cells)}
  }
}
\examples{
data(anderson)
## maybe str(anderson)
}
\keyword{datasets}



