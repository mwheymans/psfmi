\name{hipstudy}
\alias{hipstudy}
\docType{data}
\title{Dataset of elderly patients with a hip fracture}
\description{
Original dataset of elderly patients with a hip fracture
}
\usage{data(hipstudy)}
\format{
  A data frame with 426 observations on the following 18 variables.
  \describe{
    \item{\code{pat_id}}{continuous: unique patient numbers}
    \item{\code{Gender}}{dichotomous: 1 = male, 0 = female}
    \item{\code{Age}}{continuous: Years}
    \item{\code{Mobility}}{categorical: 1 = No tools, 2 = Stick / walker, 3 = Wheelchair / bed}
    \item{\code{Dementia}}{dichotomous: 2=yes, 1=no}
    \item{\code{Home}}{categorical: 1 = Independent, 2 = Elderly house, 3 = Nursering}
    \item{\code{Comorbidity}}{continuous: Number of Co_morbidities (0-4)}
    \item{\code{ASA}}{continuous: ASA score (1-4)}
    \item{\code{Hemoglobine}}{continuous: Hemoglobine pre-operative}
    \item{\code{Leucocytes}}{continuous: Leucocytes preoperative}
    \item{\code{Thrombocytes}}{continuous: Thrombocytes preoperative}
    \item{\code{CRP}}{continuous: C-reactive protein (CRP) preoperative}
    \item{\code{Creatinine}}{continuous: Creatinine preoperative}
    \item{\code{Urea}}{continuous: Urea preoperative}
    \item{\code{Albumine}}{continuous: Albumin preoperative}
    \item{\code{Fracture}}{dichotomous: 1 = per or subtrochanter fracture, 0 = collum fracture}
    \item{\code{Delay}}{continuous: time till operation in days}
    \item{\code{Mortality}}{dichotomous: 1 = yes, 0 = no}
  }
}
\examples{
data(hipstudy)
## maybe str(hipstudy)
}
\keyword{datasets}
