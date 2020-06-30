#' Function to calulate mean auc values
#' 
#'@author Martijn Heymans, 2020
#'@keywords internal
#'
#'@export
mean_auc_log <- function(auc, backtransform=TRUE){
  
  if(length(auc) ==1)
    stop("\n", "AUC contains only one value, include more AUC values", "\n")
  
  mean_auc <-
    mean(log(auc/(1-auc)))
  
  if(backtransform)
    mean_auc <-
    exp(mean_auc) /(1 + exp(mean_auc))
  
  return(mean_auc)
}
