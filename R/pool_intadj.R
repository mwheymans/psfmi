#' Provides pooled adjusted intercept after shrinkage of pooled coefficients 
#'  in multiply imputed datasets
#'
#' \code{pool_intadj} Provides pooled adjusted intercept after shrinkage of the pooled coefficients 
#'  in multiply imputed datasets for models selected with the \code{psfmi_lr} function and 
#'  internally validated with the \code{psfmi_perform} function.
#'
#' @param pobj An object of class \code{smodsmi} (selected models in multiply imputed datasets), 
#'   produced by a previous call to \code{psfmi_lr}.  
#' @param shrinkage_factor A numerical scalar. Shrinkage factor value as a result of internal validation
#'   with the \code{psfmi_perform} function.
#'  
#' @details The function provides the pooled adjusted intercept after shrinkage of pooled
#'   regression coefficients in multiply imputed datasets. The function is only available 
#'   for logistic regression models without random effects.     
#'                 
#'@return A \code{pool_intadj} object from which the following objects can be extracted: \code{int_adj}, 
#'   the adjusted intercept value, \code{coef_shrink_pooled}, the pooled regression coefficients 
#'   after shrinkage, \code{coef_orig_pooled}, the (original) pooled regression coefficients before
#'   shrinkage and \code{nimp}, the number of imputed datasets.

#' @references F. Harrell. Regression Modeling Strategies. With Applications to
#'   Linear Models, Logistic and Ordinal Regression, and Survival Analysis (2nd edition). Springer,
#'   New York, NY, 2015.
#' @references EW. Steyerberg (2019). Clinical Prediction MOdels. A Practical Approach 
#'  to Development, Validation, and Updating (2nd edition). Springer Nature Switzerland AG.
#'  
#' @references http://missingdatasolutions.rbind.io/
#' 
#' @examples
#'  res_psfmi <- psfmi_lr(data=lbpmilr, nimp=5, impvar="Impnr", Outcome="Chronic",
#'            predictors=c("Gender", "Pain","Tampascale","Smoking","Function", 
#'            "Radiation", "Age"), p.crit = 1, method="D1")
#'  res_psfmi$RR_Model
#'
#'  set.seed(100)
#'  res_val <- psfmi_perform(res_psfmi, method = "MI_boot", nboot=10, 
#'    int_val = TRUE, p.crit=1, cal.plot=FALSE, plot.indiv=FALSE)
#'  res_val$intval
#'
#'  res <- pool_intadj(res_psfmi, shrinkage_factor = 0.9774058)
#'  res$int_adj
#'  res$coef_shrink_pooled
#'   
#' @export   
pool_intadj <- function(pobj, shrinkage_factor){
  if(class(pobj)!="smodsmi")
    stop("\n", "Object should be of type smodsmi", "\n")
  if(pobj$model_type=="survival")
    stop("\n", "Pooling of intercepts only available for models of type binomial", "\n")
  if(!is.null(pobj$random.eff))
    stop("\n", "Function only available for regression models without random effects", "\n")
  if(is.null(shrinkage_factor))
    stop("\n", "Shrinkage factor not defined", "\n")
  
  data <- filter(pobj$data, get(pobj$impvar) %in% c(1:pobj$nimp))
  split_orig_imp <- data %>%
    group_by(data[,pobj$impvar])
  split_orig_data <- group_split(split_orig_imp)
  
  Y <- c(paste(pobj$Outcome, paste("~")))
  fm <- as.formula(paste(Y, paste(pobj$predictors_final, collapse = "+")))
  
  int_adj <- map(split_orig_data, function(x) {
    fit_orig <- glm(fm, data = x, family = binomial)
    coef_orig <- coef(fit_orig)[-1]
    coef_shrink <- coef_orig * shrinkage_factor
    lp_shrink <- model.matrix(fit_orig)[, -1] %*% coef_shrink
    fit_shrink <- glm(fit_orig$y ~ offset(lp_shrink), family = binomial)
    int_new <- coef(fit_shrink)
    res <- list(int_new, coef_orig, coef_shrink, names(coef_orig))
    return(res)
  })
  
  coef_orig_pooled <- colMeans(do.call("rbind",
      map(int_adj, function(x) pluck(x, 2))))
  coef_shrink_pooled <- colMeans(do.call("rbind",
      map(int_adj, function(x) pluck(x, 3))))
  int_adj <- mean(unlist(data.frame(do.call("rbind",
      map(int_adj, function(x) pluck(x, 1))))))
  
  resobj <- list(int_adj = int_adj, coef_shrink_pooled = coef_shrink_pooled,
                 coef_orig_pooled = coef_orig_pooled, nimp=pobj$nimp)
  return(resobj)
}