#' Calculates the Hosmer and Lemeshow goodness of fit test. 
#'
#' \code{hoslem_test} the Hosmer and Lemeshow goodness of fit test. 
#'
#' @param y a vector of observations (0/1).  
#' @param yhat a vector of predicted probabilities.
#' @param g Number of groups tested. Default is 10. Can not be < 3.   
#'                                               
#' @return The Chi-squared test statistic, the p-value, the observed and 
#'  expected frequencies.
#' 
#' @seealso \code{\link{pool_performance}}
#' 
#' @references Kleinman K and Horton NJ. (2014). SAS and R: Data Management, 
#'  Statistical Analysis, and Graphics. 2nd Edition. Chapman & Hall/CRC. 
#' 
#' @author Martijn Heymans, 2021
#' 
#' @examples
#'   fit <- glm(Mortality ~ Dementia + factor(Mobility) + ASA + 
#'    Gender + Age, data=hipstudy, family=binomial) 
#'    pred <- predict(fit, type = "response")
#'   
#'   hoslem_test(fit$y, pred)
#'  
#' @export 
hoslem_test <- function(y, yhat, g=10) {
  if(g <=3) stop("\n", "Number of groups must be > 3", "\n")
  cutyhat <- cut(yhat,
                breaks = quantile(yhat, probs=seq(0,
                                                  1, 1/g)), include.lowest=TRUE)
  obs <- xtabs(cbind(1 - y, y) ~ cutyhat)
  expect <- xtabs(cbind(1 - yhat, yhat) ~ cutyhat)
  chisq <- sum((obs - expect)^2/expect)
  P <- 1 - pchisq(chisq, g - 2)
  return(list(chisq=chisq, p.value=P, observed=obs, expected=expect))
}
