#' Pooling performance measures over multiply imputed datasets
#'
#' \code{pool_performance_internal} Pooling performance measures
#'
#' @param data Data frame with stacked multiple imputed datasets.
#'   The original dataset that contains missing values must be excluded 
#'   from the dataset.
#' @param formula A formula object to specify the model as normally used by glm.
#' @param nimp A numerical scalar. Number of imputed datasets. Default is 5.
#' @param impvar A character vector. Name of the variable that distinguishes 
#'   the imputed datasets.
#'
#' @keywords internal 
#'  
#' @export
pool_performance_internal <- function(data, 
                             formula, 
                             nimp, 
                             impvar)
  {
  
  coef_f <- se_f <- pred.group <- obs.group <- list()
  
  if(is_empty(formula))
    stop("\n", "Model not specified in formula object")
  
  fm <-
    formula
  
  perf_stats <-
    matrix(NA, nimp, 4)
  
  lp_mi <- matrix(NA, nrow(data[data[impvar] == 1, ]), nimp)
  
  # Pool performance measures over imputed datasets
  for (i in 1:nimp) {
    
    data_compl <- data[data[impvar] == i, ]
    f <-
      glm(fm, data = data_compl, family = binomial)
    pred <-
      predict(f, type = "response")
    
    lp_mi[, i] <- predict(f)
    
    coef_f[[i]] <-
      coef(f)
    se_f[[i]] <-
      summary(f)[[12]][, 2]
    
    # ROC/AUC
    roc.i <-
      pROC::roc(f$y, pred, quiet = TRUE)$auc
    se.roc.i <-
      sqrt(pROC::var(roc.i))
    
    # Nagelkerke R-squared
    rsq.i <-
      rsq_nagel(f)
    sc_brier.i <-
      scaled_brier(f$y, pred)
    
    perf_stats[i, ] <-
      c(roc.i, se.roc.i, rsq.i, sc_brier.i)
  }
  # End pooling performance measures in multiply imputed data

  # ROC/AUC
  roc_res <-
    pool_auc(perf_stats[, 1], perf_stats[, 2],
             nimp = nimp, log_auc = TRUE)
  
  # Pooling R square
  # Fisher z Transformation
  rsq.n <-
    tanh(mean(atanh(perf_stats[, 3])))
  
  # Colmeans of predictors in multiply imputed datasets
  coef_pooled <-
    colMeans(do.call("rbind", coef_f))
  
  sc_brier_pool <-
    mean(perf_stats[, 4])
  
  # Pooling Performance measures and coefficients
  pobjperform <- list(ROC_pooled=roc_res, coef_pooled=coef_pooled,
                      R2_pooled=rsq.n, Brier_Scaled_pooled = sc_brier_pool,
                      nimp=nimp)
  # Pooled info in each bootstrap sample
  pobjperform
}