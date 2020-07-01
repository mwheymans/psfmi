#' Calculated the scaled Brier score
#'
#' @param obs Observed outcomes.  
#' @param pred Predicted outcomes in the form of probabilities.
#'                 
#' @return The value for the scaled Brier score.
#' @seealso \code{\link{psfmi_perform}}, \code{\link{pool_performance}}
#' @author Martijn Heymans, 2020
#'   
#' @export   
scaled_brier <- function(obs, pred) {
  1 - (mean((obs - pred)^2) / (mean(obs) * (1 - mean(obs))))
}
