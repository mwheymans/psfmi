#' Function to apply RR   
#' 
#'@author Martijn Heymans, 2021
#'@keywords internal
#'
#'@export
pool_RR <- function(est, se, conf.level=0.95, n, k){
  m <- length(est)
  mean_est <- mean(est)
  var_w <-
    mean(se^2) # within variance
  var_b <-
    var(est) # between variance
  var_T <-
    var_w + (1 + (1/m)) * var_b # total variance
  se_total <-
    sqrt(var_T)
  r <- (1 + 1 / m) * (var_b / var_w)
  v_old <- (m - 1) * (1 + (1/r))^2
  lambda <- (var_b + (var_b/m))/var_T
  v_obs <- (((n-k) + 1) / ((n-k) + 3)) * (n-k) * (1-lambda)
  v_adj <- (v_old * v_obs) / (v_old + v_obs)
  alpha <- 1 - (1 - conf.level)/2
  t_stats <- mean_est/se_total
  p_val <-
    2*pt(-abs(t_stats),df=v_adj)
  t <- qt(alpha, v_adj)
  ci_upper <-
    mean_est + (t*se_total)
  ci_lower <-
    mean_est - (t*se_total)
  res <-
    round(c(mean_est, ci_lower, ci_upper, p_val), 7)
  names(res) <- c("Estimate", "95% CI L", "95% CI U", "P-val")
  return(res)
}