\name{men}
\alias{men}
\docType{data}
\title{Data of 613 patients with meningitis}
\description{
Data of 613 patients with meningitis
}
\usage{data(men)}
\format{
  A data frame with 420 observations on the following 10 variables.
  \describe{
    \item{\code{pt_id}}{continuous}
    \item{\code{sex}}{dichotomous: 0=male, 1=female}
    \item{\code{predisp}}{dichotomous: 0=no, 1=yes}
    \item{\code{mensepsi}}{categorical: disease characteristics at admission, 1=menigitis, 2=sepsis, 3=other}
    \item{\code{coma}}{dichotomous: coma at admission, 0=no, 1=coma}
    \item{\code{diastol}}{continuous: diastolic blood pressure at admission}
    \item{\code{course}}{dichotomous: disease course, 0=alive, 1=deceased}
  }
}
\examples{
data(men)
## maybe str(men)
}
\keyword{datasets}




