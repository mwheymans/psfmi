#' Function to apply RR to pool difference of NRI and AUC values  
#' 
#'@author Martijn Heymans, 2020
#'@keywords internal
#'
#'@export
RR_diff_prop <- function(est, se){
  mean_est <- mean(est)
  w_diff <-
    mean(se^2) # within variance
  b_diff <-
    var(est) # between variance
  tv_diff <-
    w_diff + (1 + (1/length(est))) * b_diff # total variance
  se_total <-
    sqrt(tv_diff)
  ci_upper <-
    mean_est + 1.96*se_total
  ci_lower <-
    mean_est - 1.96*se_total
  res <-
    c(mean_est, se_total, ci_lower, ci_upper)
  return(res)
}
