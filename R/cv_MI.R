#' Cross-validation in Multiply Imputed datasets
#'
#' \code{cv_MI} Cross-validation by applying multiple single imputation runs in train
#'  and test folds. Called by function \code{psfmi_perform}.
#'
#' @param pobj An object of class \code{pmods} (pooled models), produced by a previous
#'  call to \code{psfmi_lr}.
#' @param data_orig dataframe of original dataset that contains missing data.
#' @param folds The number of folds, default is 3.
#' @param nimp_cv Numerical scalar. Number of (multiple) imputation runs.
#' @param BW If TRUE backward selection is conducted within cross-validation. Default is FALSE.
#' @param p.crit A numerical scalar. P-value selection criterium used for backward during
#'  cross-validation. When set at 1, pooling and internal validation is done without
#'  backward selection.
#' @param miceImp Wrapper function around the \code{mice} function.
#' @param ...  Arguments as predictorMatrix, seed, maxit, etc that can be adjusted for
#'  the \code{mice} function.
#'
#' @seealso \code{\link{psfmi_perform}}
#' @author Martijn Heymans, 2020
#' @keywords internal 
#' 
#' @export
cv_MI <- function(pobj, data_orig, folds, nimp_cv, BW, p.crit, miceImp, ...)
{

  call <- match.call()

  auc_train <- auc_train_cv <- auc_test_cv <-  auc_test_se_cv <-
    coef_test <- coef_test_cv <- stats_cv_mi <- list()

  stats_cv <- matrix(NA, folds, 4)

  for(i in 1:nimp_cv){

    message("\n", "Imp run ", i, "\n")
    # Create train and test datasets
    # Stratified on outcome
    # Extract row id's in first imputed dataset
    # to apply in each imputed dataset
    idfold <-
      map(vfold_cv(data_orig, v=folds, strata = unlist(data_orig[pobj$Outcome]))$splits,
          function(x) id_test <- as.integer(row.names(x[[1]]))[-x[[2]]])

    Pred <-
      matrix(NA, nrow = nrow(data_orig), ncol = 2)

    for (f in 1:folds) {

      message("\n", "fold ", f, "\n")
      # Make copy of original dataset
      datanew <-
        data_orig
      # Outcome in test data is set to missing
      datanew[idfold[[f]], pobj$Outcome] <-
        NA

      # Apply MI (once) in the whole dataset
      imp_data <-
        miceImp(datanew, m=1, ...)
      data_compl <-
        complete(imp_data, 1)

      # Fit a model on the train data only
      Y <-
        c(paste(pobj$Outcome, paste("~")))
      fm_train <-
        as.formula(paste(Y, paste(pobj$predictors_final, collapse = "+")))

      fm_train_temp <-
        fm_train
      # if BW = TRUE
      if(BW==TRUE){
        pobj_bw <-
          bw_single(formula = fm_train_temp, data =  data_compl[-idfold[[f]], ],
                    p.crit = p.crit, keep.predictors = pobj$keep.predictors,
                    model_type="binomial")

        if(is_empty(pobj_bw$predictors_final))
          pobj_bw$predictors_final <- 1
        fm_train <-
          as.formula(paste(Y, paste(pobj_bw$predictors_final, collapse = "+")))
      }
      fit_train <-
        glm(fm_train, family = binomial, data = data_compl[-idfold[[f]], ])

      pr_train <-
        predict(fit_train, data = data_compl[-idfold[[f]], ], type="response")
      auc_train[[f]] <-
        pROC::roc(fit_train$y, pr_train, quiet = TRUE)$auc

      sc_brier_train <-
        scaled_brier(fit_train$y, pr_train)

      # Nagelkerke R-squared
      rsq_train <-
        rsq_nagel(fit_train)

      # Test the model in the test data
      Y_test <-
        unlist(data_orig[idfold[[f]], pobj$Outcome])
      pr_test <-
        predict(fit_train, newdata = data_compl[idfold[[f]], ], type="response")
      lp_test <-
        predict(fit_train, newdata = data_compl[idfold[[f]], ])
      coef_test[[f]] <-
        coef(glm(Y_test ~ lp_test, family=binomial))

      if(any(is.na(coef_test[[f]])))
        coef_test[[f]][2] <- replace_na(coef_test[[f]][2], 1)

      fit_test <-
        glm(Y_test ~ lp_test, family = binomial)
      # brier score
      sc_brier_test <-
        scaled_brier(Y_test, pr_test)
      Pred[idfold[[f]], ] <-
        c(pr_test, Y_test)

      # Nagelkerke R-squared
      rsq_test <-
        rsq_nagel(fit_test)

      stats_cv[f, ] <-
        c(sc_brier_train, sc_brier_test,
          rsq_train, rsq_test)
    }

    stats_cv_mi[[i]] <-
      stats_cv

    # AUC from package cvAUC
    cvAUC <-
      ci.cvAUC(predictions=Pred[, 1], labels=Pred[, 2],
               folds=unlist(idfold), confidence=0.95)

    auc_test_cv[[i]] <-
      cvAUC$cvAUC
    auc_test_se_cv[[i]] <-
      cvAUC$se

    auc_train_cv[[i]] <-
      auc_train

    coef_test_cv[[i]] <-
      colMeans(do.call("rbind", coef_test))
  }

  auc_train_pooled <-
    mean_auc_log(unlist(auc_train_cv))
  auc_test_pooled <-
    pool_auc(est_auc = auc_test_cv, est_se = auc_test_se_cv,
             nimp = nimp_cv, log_auc = TRUE)

  pool_stats_cv <-
    matrix(c(auc_train_pooled, auc_test_pooled[2],
             colMeans(do.call("rbind", stats_cv_mi), na.rm = TRUE)), 3,2, byrow = TRUE)
  dimnames(pool_stats_cv) <-
    list(c("AUC", "Scaled Brier", "R2"), c("Train", "Test"))

  coef_pool <-
    colMeans(do.call("rbind", coef_test_cv))

  objcv <- list(pool_stats=pool_stats_cv, LP_val=coef_pool, auc_test=auc_test_pooled)
  return(objcv)
}