#' Calculates the pooled Area Under the Curve in Multiply Imputed datasets 
#'
#' \code{pool_auc} Calculates the pooled AUC and 95% Confidence interval
#'  by using Rubin's Rules. The AUC values are log transformed before pooling. 
#'
#' @param est_auc A list of AUC values estimated in Multiply Imputed datasets.  
#' @param est_se A list of standard errors of AUC values estimated 
#'  in Multiply Imputed datasets.
#' @param nimp A numerical scalar. Number of imputed datasets. Default is 5.   
#' @param log_auc If TRUE natural logarithmic transformation is applied before
#'  pooling and finally back transformed. If FALSE the raw values are pooled.                  
#'                                               
#' @return The pooled AUC value and the 95% CI.
#' 
#' @seealso \code{\link{psfmi_perform}}, \code{\link{pool_performance}}
#' @author Martijn Heymans, 2020
#'  
#' @export   
pool_auc <- function(est_auc, est_se, nimp = 5, log_auc=TRUE){

  RR_se <- function(est, se, nimp){
    m <- nimp
    w_auc <-
      mean(se^2) # within variance
    b_auc <-
      var(est) # between variance
    tv_auc <-
      w_auc + (1 + (1/m)) * b_auc # total variance
    se_total <-
      sqrt(tv_auc)
    r <- (1 + 1 / m) * (b_auc / w_auc)
    v <- (m - 1) * (1 + (1/r))^2
    t <- qt(0.975, v)
    res <- c(se_total, t)
    return(res)
  }
  
  est_auc <-
    unlist(est_auc)
  est_auc_se <-
    unlist(est_se)
  if(length(est_auc) != nimp)
    stop("Include auc value for each imputed dataset")
  
  if(log_auc){
    est_auc_log <-
      log(est_auc/(1-est_auc))
    est_auc_se_log <-
      est_auc_se / (est_auc * (1-est_auc))
    
    se_total <-
      RR_se(est_auc_log, est_auc_se_log, nimp=nimp) # pooled se
    
    # Backtransform
    inv.auc <- exp(mean(est_auc_log)) /
      (1 + exp(mean(est_auc_log)))
    inv.auc.u <- exp(mean(est_auc_log) + (se_total[2]*se_total[1])) /
      (1 + exp(mean(est_auc_log) + (se_total[2]*se_total[1])))
    inv.auc.l <- exp(mean(est_auc_log) - (se_total[2]*se_total[1])) /
      (1 + exp(mean(est_auc_log) - (se_total[2]*se_total[1])))
    auc_res <- round(matrix(c(inv.auc.l, inv.auc, inv.auc.u),
                            1, 3, byrow = T), 4)
    dimnames(auc_res) <- list(c("AUC (logit)"),
                              c("95% Low", "AUC", "95% Up"))
  } else {
    mean_auc <-
      mean(est_auc)
    se_total <-
      RR_se(est_auc, est_auc_se, nimp=nimp)
    auc_u <-
      mean(est_auc) + (se_total[2]*se_total[1])
    if(auc_u > 1) auc_u <- 1.00
    auc_l <- mean(est_auc) - (se_total[2]*se_total[1])
    if(auc_l < 0) auc_l <- 0.00
    auc_res <-
      round(matrix(c(auc_l, mean_auc, auc_u),
                   1, 3, byrow = T), 4)
    dimnames(auc_res) <-
      list(c("AUC"), c("95% Low", "AUC", "95% Up"))
  }
  return(auc_res)
}