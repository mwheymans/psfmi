\name{chlrform}
\alias{chlrform}
\docType{data}
\title{Data about concentration of ß2-microglobuline in urine as indicator for possible damage to the kidney}
\description{
Data about concentration of ß2-microglobuline in urine as indicator for possible damage to the kidney
}
\usage{data(chlrform)}
\format{
  A data frame with 348 observations on the following 5 variables.
  \describe{
    \item{\code{pt_id}}{continuous}
    \item{\code{sport}}{categorical: degree of exposure, 0 = football player (no exposition), 1 = outdoor         swimmer (moderate exposition) and 2 = indoor swimmer (high exposition)}
    \item{\code{gammagt}}{continuous: liver damage: gamma-gt in u/l}
    \item{\code{b2}}{continuous: beta2 microglobuline in mg per mol}
    \item{\code{age}}{continuous: age in years}
  }
}
\examples{
data(chlrform)
## maybe str(chlrform)
}
\keyword{datasets}

