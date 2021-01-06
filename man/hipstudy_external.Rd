\name{hipstudy_external}
\alias{hipstudy_external}
\docType{data}
\title{External Dataset of elderly patients with a hip fracture}
\description{
External dataset of elderly patients with a hip fracture
}
\usage{data(hipstudy_external)}
\format{
  A data frame with 381 observations on the following 17 variables.
  \describe{
    \item{\code{Gender}}{dichotomous: 1 = male, 0 = female}
    \item{\code{Age}}{continuous: Years}
    \item{\code{Mobility}}{categorical: 1 = No tools, 2 = Stick / walker, 3 = Wheelchair / bed}
    \item{\code{Dementia}}{dichotomous: 2=yes, 1=no}
    \item{\code{Home}}{categorical: 1 = Independent, 2 = Elderly house, 3 = Nursering}
    \item{\code{Comorbidity}}{continuous: Number of Co-morbidities}
    \item{\code{ASA}}{continuous: ASA score}
    \item{\code{Hemoglobine}}{continuous: Hemoglobine preoperative}
    \item{\code{Leucocytes}}{continuous: Leucocytes preoperative}
    \item{\code{Thrombocytes}}{continuous: Thrombocytes preoperative}
    \item{\code{CRP}}{continuous: Creactive protein (CRP) preoperative}
    \item{\code{Creatinine}}{continuous: Creatinine preoperative}
    \item{\code{Urea}}{continuous: Urea preoperative}
    \item{\code{Albumine}}{continuous: Albumin preoperative}
    \item{\code{Fracture}}{dichotomous: 1 = per or subtrochanter fracture, 0 = collum fracture}
    \item{\code{Delay}}{continuous: time till operation in days}
    \item{\code{Mortality}}{dichotomous: 1 = yes, 0 = no}
  }
}
\examples{
data(hipstudy_external)
## maybe str(hipstudy_external)
}
\keyword{datasets}
