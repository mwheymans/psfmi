#' Function to apply RR to pool difference of NRI and AUC values  
#' 
#'@author Martijn Heymans, 2020
#'@keywords internal
#'
#'@export
RR_diff_prop <- function(est, se){
  m <- length(est)
  mean_est <- mean(est)
  w_diff <-
    mean(se^2) # within variance
  b_diff <-
    var(est) # between variance
  tv_diff <-
    w_diff + (1 + (1/m)) * b_diff # total variance
  se_total <-
    sqrt(tv_diff)
  r <- (1 + 1 / m) * (b_diff / w_diff)
  v <- (m - 1) * (1 + (1/r))^2
  t <- qt(0.975, v)
  ci_upper <-
    mean_est + t*se_total
  ci_lower <-
    mean_est - t*se_total
  res <-
    c(mean_est, se_total, ci_lower, ci_upper)
  return(res)
}
