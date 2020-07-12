#' Example dataset of Low Back Pain Patients for external validation
#'
#' Five multiply imputed datasets
#'
#' @format A data frame with 400 rows and 17 variables.
#' 
#' \describe{
#'    \item{\code{Impnr}}{a numeric vector}
#'    \item{\code{ID}}{a numeric vector}
#'    \item{\code{Chronic}}{dichotomous}
#'    \item{\code{Gender}}{dichotomous}
#'    \item{\code{Carrying}}{categorical}
#'    \item{\code{Pain}}{continuous}
#'    \item{\code{Tampascale}}{continuous}
#'    \item{\code{Function}}{continuous}
#'    \item{\code{Radiation}}{dichotomous}
#'    \item{\code{Age}}{continuous}
#'    \item{\code{Smoking}}{dichotomous}
#'    \item{\code{Satisfaction}}{categorical}
#'    \item{\code{JobControl}}{continuous}
#'    \item{\code{JobDemands}}{continuous}
#'    \item{\code{SocialSupport}}{continuous}
#'    \item{\code{Duration}}{continuous}
#'    \item{\code{BMI}}{continuous}
#'  }
#'  
#' @examples
#'  data(lbpmi_extval)
#'  ## maybe str(lbpmi_extval)\
#'  
#' @keywords dataset
#' 
"lbpmi_extval"