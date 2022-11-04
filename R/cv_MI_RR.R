#' Cross-validation in Multiply Imputed datasets
#'
#' \code{cv_MI_RR} Cross-validation by applying multiply imputed pooled models in train
#'  and test folds. Called by function \code{psfmi_perform}.
#'
#' @param pobj An object of class \code{pmods} (pooled models), produced by a previous
#'  call to \code{psfmi_lr}.
#' @param data_orig dataframe of original dataset that contains missing data.
#' @param folds The number of folds, default is 3.
#' @param nimp_mice Numerical scalar. Number of multiple imputation runs.
#' @param p.crit A numerical scalar. P-value selection criterium used for backward during
#'  cross-validation. When set at 1, pooling and internal validation is done without
#'  backward selection.
#' @param BW If TRUE backward selection is conducted within cross-validation. Default is FALSE.
#' @param miceImp Wrapper function around the \code{mice} function.
#' @param ...  Arguments as predictorMatrix, seed, maxit, etc that can be adjusted for
#'  the \code{mice} function.
#'
#' @seealso \code{\link{psfmi_perform}}
#' @author Martijn Heymans, 2020
#' @keywords internal 
#' 
#' @export
cv_MI_RR <- function(pobj, data_orig, folds, nimp_mice, p.crit, BW, miceImp, ...)
{

  # Single fold definition
  idfold <- map(vfold_cv(data_orig, v=folds,
                         strata = pobj$Outcome)$splits,
                function(x)
                  id_test <- as.integer(row.names(x[[1]]))[-x[[2]]])

  Y <- c(paste(pobj$Outcome, paste("~")))
  fm_train <- as.formula(paste(Y, paste(pobj$predictors_final, collapse = "+")))

  fm_train_temp <- fm_train

  modeldim <-
    dim(model.matrix(fm_train , data = data_orig))[2]

  coefs_RR <-
    matrix(NA, nrow = folds, ncol = modeldim)

  coef_m_train <-
    matrix(NA, nrow = nimp_mice, ncol = modeldim)

  Xb <-
    Pred <-
    matrix(NA, nrow = nrow(data_orig), ncol = nimp_mice)

  auc_train_mi <- se_auc_train_mi <- sc_brier_train_mi <-
    rsq_train_mi <- rsq_train_RR <- sc_brier_train_RR <-
    auc_train_RR <- lp_test_RR <- fit_test_RR_rsq <-
    sc_brier_test_RR <- auc_test_f <- auc_test_RR <- rsq_test_RR <- list()

  for (f in 1:folds) {

    message("\n", "fold ", f, "\n")

    datanew <-
      data_orig
    # Outcome in test data is set to missing
    datanew[idfold[[f]], pobj$Outcome] <-
      NA

    X_m_test <-
      array(NA, dim = c(length(idfold[[f]]), modeldim, nimp_mice))

    # apply mice
    imp_data <-
      miceImp(datanew, m=nimp_mice, ...)

    imp_data_compl <- complete(imp_data, action = "long", include = FALSE)

    # if BW = TRUE
    if(BW==TRUE){
      split_data <- group_split(imp_data_compl %>%
                                  group_by(imp_data_compl[, ".imp"]))
      # select in each imputed dataset the same bootstrap id's
      imp_data_compl <-
        do.call("rbind", lapply(split_data, function(x) x[-idfold[[f]], ]))

      pobj_bw <- psfmi_lr(formula = fm_train_temp, data =  imp_data_compl, nimp=nimp_mice, impvar = ".imp",
                          p.crit = p.crit, keep.predictors = pobj$keep.predictors,
                          method = pobj$method, direction = "BW")

      if(is_empty(pobj_bw$predictors_final))
        pobj_bw$predictors_final <- 1

      fm_train <-
        as.formula(paste(Y, paste(pobj_bw$predictors_final, collapse = "+")))

      modeldim <-
        dim(model.matrix(fm_train , data = data_orig))[2]

      coefs_RR <- matrix(NA, nrow = folds, ncol = modeldim)

      coef_m_train <-
        matrix(NA, nrow = nimp_mice, ncol = modeldim)

      X_m_test <-
        array(NA, dim = c(length(idfold[[f]]), modeldim, nimp_mice))
    }

    for (m in 1:nimp_mice) {
      data_compl <-
        complete(imp_data, m)  # Select the completed data

      fit_train <-
        glm(fm_train, family = binomial, data = data_compl[-idfold[[f]], ])

      # Save performance of models in training sets
      pr_train <-
        predict(fit_train, type="response")
      auc_train_mi[[m]] <-
        pROC::roc(fit_train$y, pr_train, quiet = TRUE )$auc
      se_auc_train_mi[[m]] <- sqrt(pROC::var(auc_train_mi[[m]]))
      sc_brier_train_mi[[m]] <-
        scaled_brier(fit_train$y, pr_train)

      rsq_train_mi[[m]] <- rsq_nagel(fit_train)

      # Set aside X of completed test data
      X_m_test[, , m]  <-
        model.matrix(fm_train,
                     data = data_compl[idfold[[f]], ],
                     drop.unused.levels = FALSE)
      # Save the train coefficients
      coef_m_train[m, ] <- fit_train$coefficients
    }

    # Pool train data performance measures from imputed datasets
    sc_brier_train_RR[[f]] <-
      mean(unlist(sc_brier_train_mi))
    auc_train_RR[[f]] <-
      pool_auc(est_auc = unlist(auc_train_mi),
               est_se = unlist(se_auc_train_mi),
               nimp = nimp_mice, log_auc = TRUE)[2]

    rsq_train <- atanh(unlist(rsq_train_mi))
    z.rsq.train <- mean(rsq_train)
    rsq_train_RR[[f]] <- tanh(z.rsq.train)

    # Pool the train coefficients
    coefs_RR[f, ] <-
      apply(coef_m_train, 2, mean)

    # Apply pooled coefficients to
    # all imputed validation data matrices
    # and save linear predictors
    Xb[idfold[[f]], ] <-
      apply((X_m_test) * ((rep(
        1, length(idfold[[f]])
      )) %*% t(coefs_RR[f,])) %o% (rep(1, nimp_mice)), c(1, 3), sum)

    fit_train$coefficients <- coefs_RR[f, ]

    if (!any(is.na(data_orig[idfold[[f]], ]))) {
      # complete record
      Pred[idfold[[f]], ] <-
        matrix(predict.glm(fit_train, data_compl[idfold[[f]], ], type = "response"), ncol =
                 1) %*% rep(1, nimp_mice)
    } else {
      for (m in 1:nimp_mice) {
        data_compl <- complete(imp_data, m)
        Pred[idfold[[f]], m] <-
          predict.glm(fit_train, data_compl[idfold[[f]], ], type = "response")
      }
    }

    # Pool test data performance measures from imputed datasets
    Obs_test_outcome <-
      unlist(data_orig[idfold[[f]], pobj$Outcome])
    Pred_test_outcome <-
      Pred[idfold[[f]], ]

    sc_brier_test_RR[[f]] <-
      mean(apply(Pred_test_outcome, 2, function(x)
        scaled_brier(Obs_test_outcome, x) ))
    # Test the model in the test data
    LP_test <-
      Xb[idfold[[f]], ]

    lp_test_RR[[f]] <-
      colMeans(t(apply(LP_test, 2, function(x)
        coef(glm(Obs_test_outcome ~ x, family=binomial)))))

    auc_test_f[[f]] <-
      apply(Pred_test_outcome, 2, function(x){
        roc_test <- pROC::roc(Obs_test_outcome, x, quiet = TRUE)$auc
        roc_test_se <- se_auc_train_mi[[m]] <- sqrt(pROC::var(roc_test))
        c(roc_test, roc_test_se)
      })
    auc_test_RR[[f]] <-
      pool_auc(est_auc = as.list(auc_test_f[[f]][1, ]),
               est_se = as.list(auc_test_f[[f]][2, ]), nimp = nimp_mice, log_auc = TRUE)[2]

    fit_test_RR_rsq[[f]] <-
      apply(LP_test, 2, function(x) {
        fit_test_rsq <- glm(Obs_test_outcome ~ x, family=binomial)
        rsq.nagel <- rsq_nagel(fit_test_rsq)
      })

    rsq_test <- atanh(unlist(fit_test_RR_rsq))
    z.rsq.test <- mean(rsq_test)

    rsq_test_RR[[f]] <- tanh(z.rsq.test)

  }

  avg_sc_brier_train <-
    mean(unlist(sc_brier_train_RR))
  avg_sc_brier_test <-
    mean(unlist(sc_brier_test_RR))

  avg_auc_train <-
    mean_auc_log(unlist(auc_train_RR))
  avg_auc_test <-
    mean_auc_log(unlist(auc_test_RR))
  avg_rsq_train <-
    mean(unlist(rsq_train_RR))
  avg_rsq_test <-
    mean(unlist(rsq_test_RR))

  avg_LP_test <-
    colMeans(do.call("rbind", lp_test_RR))
  names(avg_LP_test) <-
    c("Intercept", "Slope")

  cv_mi_stats <-
    data.frame(matrix(c(avg_auc_train, avg_auc_test,
                        avg_sc_brier_train, avg_sc_brier_test,
                        avg_rsq_train, avg_rsq_test), 3, 2, byrow = TRUE))
  row.names(cv_mi_stats) <-
    c("AUC", "Brier scaled", "Rsq")
  names(cv_mi_stats) <-
    c("Train", "Test")

  objcv <- list(stats=cv_mi_stats, slope=avg_LP_test)
  objcv
}