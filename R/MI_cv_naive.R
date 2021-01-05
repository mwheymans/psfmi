#' Naive method for Cross-validation in Multiply Imputed datasets
#'
#' \code{MI_cv_naive} Cross-validation by applying multiply imputed pooled models in train
#'  and test folds. Called by function \code{psfmi_perform}.
#'
#' @param pobj An object of class \code{pmods} (pooled models), produced by a previous
#'  call to \code{psfmi_lr}.
#' @param folds The number of folds, default is 3.
#' @param p.crit A numerical scalar. P-value selection criterium used for backward during
#'  cross-validation. When set at 1, pooling and internal validation is done without
#'  backward selection.
#' @param BW If TRUE backward selection is conducted within cross-validation. Default is FALSE.
#' @param cv_naive_appt Default is TRUE for showing the cross-validation apparent (train) and
#' test results. Set to FALSE to only give test results.
#' 
#' @seealso \code{\link{psfmi_perform}}
#' @author Martijn Heymans, 2020
#' @keywords internal 
#' 
#' @export
MI_cv_naive <- function(pobj, folds = 3, p.crit = 1, BW=FALSE, cv_naive_appt=TRUE)
{
  # Create train and test datasets
  # Stratified on outcome
  # Extract row id's in first imputed datseet
  # to apply in each imputed dataset
  idfold <- map(vfold_cv(pobj$data[pobj$data[pobj$impvar] == 1, ],
                         v=folds, strata = pobj$Outcome)$splits,
                function(x) {
                  id_train <- x[[2]]
                  id_test <- as.integer(row.names(x[[1]]))[-id_train]
                })

  # Apply model in train and test set
  # and in each imputed dataset

  test_cv <- auc.i_train <- auc_se.i_train <-
    rsq.i_train <- sc_brier.i_train <- cal_coef.i <-
    auc.i_test <- auc_se.i_test <- rsq.i_test <-
    sc_brier.i_test <- cal_coef.i <- list()
  for(i in 1:pobj$nimp)
  {
    message("\n", "Imputation ", i)
    dat_imp <-
      pobj$data[pobj$data[pobj$impvar] == i, ]

    cv_perform <- lapply(idfold,
                         function(x) {

                           # Apply model in train data
                           train_data <- dat_imp[-x, ]
                           Y <- c(paste(pobj$Outcome, paste("~")))
                           fm_train <- as.formula(paste(Y, paste(pobj$predictors_final, collapse = "+")))

                           # if BW = TRUE
                           if(BW==TRUE){
                             pobj_bw <- bw_single(formula = fm_train, data =  train_data,
                                                  p.crit = p.crit, keep.predictors = pobj$keep.predictors)

                             if(is_empty(pobj_bw$predictors_final))
                               pobj_bw$predictors_final <- 1
                             fm_train <-
                               as.formula(paste(Y, paste(pobj_bw$predictors_final, collapse = "+")))
                           }

                           fit_train <- glm(fm_train, data=train_data, family=binomial)

                           # Apply model in train data
                          pr_train <-
                             predict(fit_train, type="response")

                           auc_train <- pROC::roc(fit_train$y, pr_train, quiet = TRUE)$auc
                           se_auc_train <- sqrt(pROC::var(auc_train))

                           sc_brier_train <-
                             scaled_brier(fit_train$y, pr_train)

                           # Nagelkerke R-squared
                           rsq_train <- rsq_nagel(fit_train)

                           # Apply model in test data
                           # Apply model in test data
                           test_data <- dat_imp[x, ]
                           pr_test <- predict(fit_train, newdata = test_data, type="response")
                           lp_test <- predict(fit_train, newdata = test_data)
                           coef_test <- coef(glm(unlist(test_data[pobj$Outcome]) ~ lp_test, family=binomial))
                           fit_test <- glm(unlist(test_data[pobj$Outcome]) ~ lp_test, family=binomial)

                           if(any(is.na(coef_test)))
                             coef_test[2] <- replace_na(coef_test[2], 1)

                           # brier score
                           sc_brier_test <- scaled_brier(fit_test$y, pr_test)

                           # Nagelkerke R-squared
                           rsq_test <- rsq_nagel(fit_test)


                           list(folds=x, pred_train=pr_train, obs_train=fit_train$y,
                                pred_test=pr_test, obs_test=fit_test$y,
                                rsq_train=rsq_train, rsq_test=rsq_test,
                                sc_brier_train=sc_brier_train,
                                sc_brier_test=sc_brier_test,
                                coef_test=coef_test, auc_train = auc_train, se_auc_train)
                         })

    id_folds_test <-
      sapply(cv_perform, function(x) x[1])
    pred_outcome_test <-
      unlist(sapply(cv_perform, function(x) x[4]))
    obs_outcome_test <-
      unlist(sapply(cv_perform, function(x) x[5]))

    # Take mean of logit AUC
    auc.i_train[[i]] <-
      mean_auc_log(sapply(cv_perform, function(x) x$auc_train))

    cvAUC_test <-
      ci.cvAUC(predictions=pred_outcome_test, labels=obs_outcome_test,
               folds=id_folds_test, confidence=0.95)
    auc.i_test[[i]] <-
      cvAUC_test$cvAUC
    auc_se.i_test[[i]] <-
      cvAUC_test$se

    rsq.i_train[[i]] <-
      mean(unlist(sapply(cv_perform,
                         function(x) x$rsq_train)), na.rm = TRUE)
    sc_brier.i_train[[i]] <-
      mean(unlist(sapply(cv_perform,
                         function(x) x$sc_brier_train)), na.rm = TRUE)

    rsq.i_test[[i]] <-
      mean(unlist(sapply(cv_perform,
                         function(x) x$rsq_test)), na.rm = TRUE)
    sc_brier.i_test[[i]] <-
      mean(unlist(sapply(cv_perform,
                         function(x) x$sc_brier_test)), na.rm = TRUE)
    cal_coef.i[[i]] <-
      colMeans(t(sapply(cv_perform,
                        function(x) x$coef_test)), na.rm = TRUE)
  }


  # Pooling R square
  # Fisher z Transformation
  z.rsq_train <-
    atanh(unlist(rsq.i_train))
  z.rsq.p_train <-
    mean(z.rsq_train)

  # inv Fisher z = pooled rsq
  pool_R2_train <-
    tanh(z.rsq.p_train)
  pool_sc_brier_train <-
    mean(unlist(sc_brier.i_train))

  # Pooling R square
  # Fisher z Transformation
  z.rsq_test <-
    atanh(unlist(rsq.i_test))
  z.rsq.p_test <-
    mean(z.rsq_test)

  # inv Fisher z = pooled rsq
  pool_R2_test <-
    tanh(z.rsq.p_test)
  pool_sc_brier_test <-
    mean(unlist(sc_brier.i_test))

  pool_coef <-
    colMeans(do.call("rbind",  cal_coef.i) )
  names(pool_coef) <-
    c("Intercept", "Slope")

  auc_pooled_train <-
    mean_auc_log(unlist(auc.i_train))
  auc_pooled_test <-
    pool_auc(est_auc = auc.i_test, est_se = auc_se.i_test, nimp = pobj$nimp, log_auc = TRUE)

  auc <-
    c(auc_pooled_train, auc_pooled_test[2])
  sc_brier <-
    c(pool_sc_brier_train, pool_sc_brier_test)
  rsq <-
    c(pool_R2_train, pool_R2_test)

  cv_stats <-
    data.frame(matrix(c(auc, sc_brier, rsq), 3, 2, byrow = TRUE))
  row.names(cv_stats) <-
    c("AUC", "Brier scaled", "R-squared")
  names(cv_stats) <-
    c("Train", "Test")

  rescv <-
    list(cv_stats = cv_stats, auc_test=auc_pooled_test, test_coef=pool_coef)

  if(cv_naive_appt){
    Y <- c(paste(pobj$Outcome, paste("~")))
    if(is_empty(pobj$predictors_final)) {
      pobj$predictors_final <- 1
      fm <- as.formula(paste(Y, paste(pobj$predictors_final, collapse = "+")))
    } else {
      fm <- as.formula(paste(Y, paste(pobj$predictors_final, collapse = "+")))
    }

    perform_mi_orig <-
      pool_performance(data=pobj$data, nimp = pobj$nimp,
                       impvar=pobj$impvar, Outcome = pobj$Outcome,
                       predictors = pobj$predictors_final, cal.plot=FALSE,
                       plot.indiv=FALSE, groups_cal = 10)
    rescv <- list(Test=rescv, Apparent=perform_mi_orig)
  }
  return(rescv)
}