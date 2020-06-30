#' Bootstrap validation in Multiply Imputed datasets
#'
#' \code{MI_boot} Bootstrapping in each (original) Multiply Imputed dataset for internal validation.
#'  Called by function \code{psfmi_perform}.
#'
#' @param pobj An object of class \code{pmods} (pooled models), produced by a previous
#'  call to \code{psfmi_lr}.
#' @param p.crit A numerical scalar. P-value selection criterium used for backward
#'  or forward selection during in the bootstrap samples. When set at 1, validation is done
#'  without variable selection.
#' @param nboot The number of bootstrap resamples, default is 10.
#' @param direction The direction of predictor selection, "BW" is for backward selection and
#'  "FW" for forward selection.
#'  
#' @seealso \code{\link{psfmi_perform}}
#' @author Martijn Heymans, 2020
#' @keywords internal 
#' 
#' @export
MI_boot <- function(pobj, p.crit, nboot, direction)
{
  call <- match.call()

  imp1 <-
    pobj$data[pobj$data[pobj$impvar] == 1, ]
  boot_data <-
    bootstraps(imp1, times = nboot)

  boot_seq <-
    as.list(1:nboot)

  opt_boot <-
    mapply(function(x, y) {

      message("\n", "Boot ", y)
      boot_id <-
        x[[2]]

      # Split original imputed data and select
      # same bootstrap cases in imputed datasets
      split_orig_imp <-
        pobj$data %>% group_by(pobj$data[,pobj$impvar])
      split_app_data <-
        group_split(split_orig_imp)
      # select in each imputed dataset the same bootstrap id's
      boot_app_data <-
        lapply(split_app_data, function(x) {
        x <- as.data.frame(x)
        app_data_select <- x[boot_id, ]
      })
      app_data <-
        do.call("rbind", boot_app_data)

      Y_boot <-
        c(paste(pobj$Outcome, paste("~")))

      fm_boot <-
        as.formula(paste(Y_boot, paste(pobj$predictors_final, collapse = "+")))

      pool_model_lr <-
        psfmi_lr(formula = fm_boot, data =  app_data, nimp=pobj$nimp, impvar = pobj$impvar,
                                p.crit = p.crit, keep.predictors = pobj$keep.predictors,
                                method = pobj$method, direction = direction)
      if(p.crit!=1){
        if(direction=="BW")
          predictors_selected <-
            ifelse(pool_model_lr$predictors_out[nrow(pool_model_lr$predictors_out), ], 0, 1)
        if(direction=="FW")
          predictors_selected <-
            pool_model_lr$predictors_in[nrow(pool_model_lr$predictors_in), ]
      } else {
        predictors_selected <-
          ifelse(pool_model_lr$predictors_out[nrow(pool_model_lr$predictors_out), ], 0, 1)
      }

      Y <- c(paste(pobj$Outcome, paste("~")))
      if(is_empty(pool_model_lr$predictors_final)) {
        pool_model_lr$predictors_final <- 1
        fm <-
          as.formula(paste(Y, paste(pool_model_lr$predictors_final, collapse = "+")))
        lp_app_pooled <- 1
      } else {
        fm <-
          as.formula(paste(Y, paste(pool_model_lr$predictors_final, collapse = "+")))
        lp_app_pooled <- pool_model_lr$RR_model[[1]][, 2]
        if(p.crit<1)
          lp_app_pooled <- pool_model_lr$RR_model_final[[1]][, 2]
      }

      # Obtain apparent pooled performance measures
      perform_app <-
        pool_performance(data=app_data, nimp = pobj$nimp,
                         impvar=pobj$impvar, Outcome = pobj$Outcome,
                         predictors = pool_model_lr$predictors_final,
                         cal.plot = FALSE)

      # Test apparent LP in (imputed) original data
      coef_slope_test <-
        list()

      perform_test <-
        matrix(NA, pobj$nimp, 4)

      for(i in 1:pobj$nimp){
        data_test <-
          pobj$data[pobj$data[pobj$impvar] == i, ]

        fit <-
          glm(fm, data = data_test, family = binomial)

        lp_test <-
          model.matrix(fit) %*% lp_app_pooled
        fit_test <-
          glm(fit$y ~ lp_test, family = binomial)
        coef_fit_test <-
          coef(fit_test)

        if(length(coef(fit))==1)
          coef_fit_test <- replace_na(coef_fit_test, 1)
        coef_slope_test[[i]] <- coef_fit_test

        test_prob <-
          c(1/(1+exp(-lp_test)))

        # ROC/AUC
        roc_test <-
          roc(fit$y, test_prob, quiet = TRUE)$auc
        roc_test_se <-
          sqrt(pROC::var(roc_test))

        # Nagelkerke R-squared
        rsq_test <-
          rsq_nagel(fit_test)

        # Brier and scaled Brier score
        sc_brier_test <-
          scaled_brier(fit$y, test_prob)

        perform_test[i, ] <-
          c(roc_test, roc_test_se,
            rsq_test, sc_brier_test)
      }

      roc_test_pool <-
        pool_auc(perform_test[, 1], perform_test[, 2], nimp = pobj$nimp, log_auc = TRUE)

      # Pooling R square
      # Fisher z Transformation
      z.rsq <- atanh(perform_test[, 3])
      z.rsq.p <- mean(z.rsq)

      # within variance
      n <- nrow(pobj$data[pobj$data[pobj$impvar] == 1, ])
      se.z.rsq <- 1/(n-3)
      # between variance
      b.rsq <- var(z.rsq)
      # total variance
      tv.rsq <- se.z.rsq + ((1 + (1/pobj$nimp)) * b.rsq)
      se.t.rsq <- sqrt(tv.rsq)
      # inv Fisher z = pooled rsq
      R2_test <- tanh(z.rsq.p)

      sc_brier_pool_test <- mean(perform_test[, 4])

      lp_test_pooled <- colMeans(do.call("rbind", coef_slope_test), na.rm = TRUE)
      # End pooling performance measures in multiply imputed data

      ROC_app <-
        perform_app$ROC_pooled[2]
      ROC_test <-
        roc_test_pool[2]
      Slope_test <-
        lp_test_pooled
      R2_app <-
        perform_app$R2_pooled
      Brier_scaled_app <-
        perform_app$Brier_Scaled_pooled
      Brier_scaled_test <-
        sc_brier_pool_test

      opt_perform <-
        list(c(ROC_app, ROC_test, R2_app, R2_test,
               Brier_scaled_app, Brier_scaled_test, Slope_test),
             predictors_selected)
      return(opt_perform)

    }, x = boot_data$splits, y=boot_seq, SIMPLIFY = FALSE)

  predictors_selected <-
    data.frame(do.call("rbind", lapply(opt_boot, function(x) x[[2]])))
  colnames(predictors_selected) <-
    pobj$predictors_initial
  row.names(predictors_selected) <-
    paste("Boot", 1:nboot)

  res_boot <-
    data.frame(do.call("rbind", lapply(opt_boot, function(x) x[[1]])))
  colnames(res_boot) <-
    c("ROC_app", "ROC_test", "R2_app", "R2_test",
      "Brier_sc_app", "Brier_sc_test",
      "intercept", "Slope")
  row.names(res_boot) <-
    paste("Boot", 1:nboot)

  roc_optimism <-
    res_boot$ROC_app - res_boot$ROC_test
  r2_optimism <-
    res_boot$R2_app - res_boot$R2_test
  sc_brier_optimism <-
    res_boot$Brier_sc_app - res_boot$Brier_sc_test

  res_boot_m <-
    colMeans(data.frame(res_boot,
                        roc_optimism, r2_optimism, sc_brier_optimism), na.rm=TRUE)

  # Perform original model in multiply imputed original data
  Y_orig <- c(paste(pobj$Outcome, paste("~")))
  if(is_empty(pobj$predictors_final)) {
    pobj$predictors_final <- 1
    fm_orig <-
      as.formula(paste(Y_orig, paste(pobj$predictors_final, collapse = "+")))
  } else {
    fm_orig <-
      as.formula(paste(Y_orig, paste(pobj$predictors_final, collapse = "+")))
  }
  perform_mi_orig <-
    pool_performance(data=pobj$data, nimp = pobj$nimp,
                     impvar=pobj$impvar, Outcome = pobj$Outcome,
                     predictors = pobj$predictors_final, cal.plot = FALSE)

  ROC_orig <-
    perform_mi_orig$ROC_pooled[2]
  R2_orig <-
    perform_mi_orig$R2_pooled
  sc_Brier_orig <-
    perform_mi_orig$Brier_Scaled_pooled

  ROC_corr <-
    ROC_orig - res_boot_m["roc_optimism"]
  ROC_val <-
    c(ROC_orig, res_boot_m["ROC_app"],
      res_boot_m["ROC_test"], res_boot_m["roc_optimism"], ROC_corr)

  R2_corr <-
    R2_orig - res_boot_m["r2_optimism"]
  R2_val <-
    c(R2_orig, res_boot_m["R2_app"],
      res_boot_m["R2_test"], res_boot_m["r2_optimism"], R2_corr)

  sc_Brier_corr <-
    sc_Brier_orig - res_boot_m["sc_brier_optimism"]
  sc_Brier_val <-
    c(sc_Brier_orig, res_boot_m["Brier_sc_app"],
      res_boot_m["Brier_sc_test"], res_boot_m["sc_brier_optimism"], sc_Brier_corr)

  Slope_corr <-
    res_boot_m["Slope.test"]
  Slope_val <-
    c(1, 1, res_boot_m["Slope"], 1-res_boot_m["Slope"], res_boot_m["Slope"])

  pobjval <-
    as.data.frame(matrix(c(ROC_val, R2_val, sc_Brier_val, Slope_val), 4, 5, byrow = T))

  colnames(pobjval) <-
    c("Orig", "Apparent", "Test", "Optimism", "Corrected")
  row.names(pobjval) <-
    c("AUC", "R2", "Brier Scaled", "Slope")
  pobjval <- list(stats_val = pobjval, intercept_test = res_boot_m[7], res_boot = res_boot)
  return(pobjval)
}
