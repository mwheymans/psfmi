\name{infarct}
\alias{infarct}
\docType{data}
\title{Data of a patient-control study regarding the relationship between MI and smoking}
\description{
Data of a patient-control study regarding the relationship between MI and smoking
}
\usage{data(infarct)}
\format{
  A data frame with 420 observations on the following 10 variables.
  \describe{
    \item{\code{ppnr}}{continuous}
    \item{\code{infarct}}{dichotomous: 1=yes, 0=no}
    \item{\code{smoking}}{dichotomous: 1=yes, 0=no}
    \item{\code{alcohol}}{categorical: 0= none, 1= 1-2 glasses a day, 2= >2 glasses a day}
    \item{\code{active}}{dichotomous: 1=active, 0=inactive}
    \item{\code{sex}}{dichotomous: 1=male, 0=female}
    \item{\code{profession}}{categorical:1=epidemiologist, 2=statistician, 3=other}
    \item{\code{bmi}}{continuous: body mass index}
    \item{\code{sys}}{continuous: systolic blood pressure}
    \item{\code{dias}}{continuous: diastolic blood pressure}
  }
}
\examples{
data(infarct)
## maybe str(infarct)
}
\keyword{datasets}
