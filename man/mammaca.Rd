\name{mammaca}
\alias{mammaca}
\docType{data}
\title{Data of a study among women with breast cancer}
\description{
Data of a study among women with breast cancer
}
\usage{data(mammaca)}
\format{
  A data frame with 1207 observations on the following 10 variables.
  \describe{
    \item{\code{id}}{continuous}
    \item{\code{time}}{continuous, Time (months)}
    \item{\code{status}}{dichotomous: 1=yes, 0=no}
    \item{\code{er}}{Estrogen Receptor Status, 1=positive, 0=negative}
    \item{\code{age}}{continuous}
    \item{\code{histgrad}}{categorical}
    \item{\code{ln_yesno}}{lymph nodes, 0=no, 1=yes}
    \item{\code{pathsd}}{dichotomous: Pathological Tumor Size}
    \item{\code{pr}}{dichotomous: Progesterone Receptor Status, 0=negative, 1=positive}
  }
}
\examples{
data(mammaca)
## maybe str(mammaca)
}
\keyword{datasets}
