#' Combines the Chi Square statistics across Multiply Imputed datasets 
#'
#' \code{pool_D2} The D2 statistic to combine the Chi square values 
#'  across Multiply Imputed datasets.
#'  
#' @param dw a vector of Chi square values obtained after multiple imputation.  
#' @param v single value for the degrees of freedom of the Chi square statistic.
#'                                               
#' @return The pooled chi square values as the D2 statistic, the p-value,
#'  the numerator, df1 and denominator, df2 degrees of freedom for the
#'  F-test.
#' 
#' @references Eekhout I, van de Wiel MA, Heymans MW. Methods for significance testing of categorical
#'   covariates in logistic regression models after multiple imputation: power and applicability
#'   analysis. BMC Med Res Methodol. 2017;17(1):129.
#' @references Van Buuren S. (2018). Flexible Imputation of Missing Data. 2nd Edition. Chapman & Hall/CRC
#'   Interdisciplinary Statistics. Boca Raton.
#'   
#' @author Martijn Heymans, 2021
#'  
#' @examples
#'   pool_D2(c(2.25, 3.95, 6.24, 5.27, 2.81), 4) 
#'  
#' @export 
pool_D2 <- function(dw, v){
  m <-
    length(dw)
  dw_sqrt <-
    sqrt(dw)
  r <-
    (1 + 1/m) * var(dw_sqrt)
  D2 <-
    (mean(dw)/v - (m + 1)/(m - 1) * r) / (1 + r)
  v2 <-
    1/v^(3/m) * ((m - 1) * (1 + 1/r)^2)
  p_value <-
    stats::pf(D2, df1 = v, df2 = v2, lower.tail = FALSE)
  res <-
    matrix(c(D2, p_value, v, v2), ncol=4)
  colnames(res) <- c("F_value", "P(>F)", "df1", "df2")
  return(res)
}
