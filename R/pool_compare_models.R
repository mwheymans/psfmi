#' Compare the fit and performance of prediction models across Multipy Imputed data 
#'
#' \code{pool_compare_model} Compares the fit and performance of prediction models 
#'  in multiply imputed data sets by using clinical important performance measures
#'
#' @param pobj An object of class \code{pmods} (pooled models), produced by a previous
#'  call to \code{psfmi_lr}.
#' @param compare.predictors Character vector with the names of the predictors that are 
#'  compared. See details.
#' @param compare.group Character vector with the names of the group of predictors that are 
#'  compared. See details.
#' @param cutoff A numerical scalar. Cutoff used for the categorical NRI value. More than one
#'  cutoff value can be used.
#' @param boot_auc If TRUE the standard error of the AUC is calculated with stratified
#'  bootstrapping. If FALSE (is default), the standard error is calculated with De Long's
#'  method.
#' @param nboot A numerical scalar. The number of bootstrap samples for the AUC standard error, used 
#'  when boot_auc is TRUE. Default is 1000.      
#'  
#' @details The fit of the models are compared by using the D3 method for pooling Likelihood ratio 
#'  statistics (method of Meng and Rubin). The pooled AIC difference is calculated according to
#'  the formula \code{AIC = D - 2*p}, where D is the pooled likelihood ratio tests of 
#'  constrained models (numerator in D3 statistic) and p is the difference in number of parameters 
#'  between the full and restricted models that are compared. The pooled AUC difference  
#'  is calculated, after the standard error is obtained in each imputed data set by method 
#'  DeLong or bootstrapping. The NRI categorical and continuous and IDI are calculated in each 
#'  imputed data set and pooled.
#'    
#'@return An object from which the following objects can be extracted: 
#'  \itemize{
#'  \item  \code{DR_stats} p-value of the D3 statistic, the D3 statistic, LRT fixed is the 
#'   likelihood Ratio test value of the constrained models.
#'  \item  \code{stats_compare} Mean of LogLik0, LogLik1, AIC0, AIC1, AIC_diff values of the 
#'   restricted (containing a 0) and full models (containing a 1).
#'  \item  \code{NRI} pooled values for the categorical and continuous Net Reclassification
#'   improvement values and the Integrated Discrimination improvement.
#'  \item  \code{AUC_stats} Pooled Area Under the Curve of restricted and full models. 
#'  \item  \code{AUC_diff} Pooled difference in AUC. 
#'  \item  \code{formula_test} regression formula of full model.
#'  \item  \code{cutoff} Cutoff value used for reclassification values.
#'  \item  \code{formula_null} regression formula of null model
#'  \item  \code{compare_predictors} Predictors used in full model.
#'  \item  \code{compare_group} group of predictors used in full model.
#' }
#'  
#' @references Eekhout I, van de Wiel MA, Heymans MW. Methods for significance testing of categorical
#'   covariates in logistic regression models after multiple imputation: power and applicability
#'   analysis. BMC Med Res Methodol. 2017;17(1):129.
#' @references Consentino F, Claeskens G. Order Selection tests with multiply imputed data
#'   Computational Statistics and Data Analysis.2010;54:2284-2295.   
#'  
#' @examples
#'  pool_lr <- psfmi_lr(data=lbpmilr, p.crit = 1, direction="FW", nimp=10, impvar="Impnr", 
#'  Outcome="Chronic", predictors=c("Radiation"), cat.predictors = ("Satisfaction"),
#'  int.predictors = NULL, spline.predictors="Tampascale", nknots=3, method="D1")
#'
#'  res_compare <- pool_compare_models(pool_lr, compare.predictors = c("Pain", "Duration", 
#'  "Function"), cutoff = 0.4)
#'  res_compare
#'
#'  
#' @export
pool_compare_models <- function(pobj,
                               compare.predictors = NULL,
                               compare.group = NULL,
                               cutoff = 0.5,
                               boot_auc=FALSE,
                               nboot=1000)
{
  
  call <- match.call()
  
  if(!inherits(pobj, "pmods"))
    stop("\n", "Object should be of type pmods", "\n")
  if(pobj$model_type=="survival")
    stop("\n", "Methods only available for models of type binomial", "\n")
  if(is_null(compare.predictors) & is_null(compare.group))
    stop("\n", "Define at leat 1 predictor or group of predictors to compare", "\n")
  if(length(compare.group)==1)
    stop("\n", "compare.group must contain > 1 predictor, otherwise use compare.predictors", "\n")
  
  P_temp_final <- pobj$predictors_final
  P_temp_final <- clean_P(P_temp_final)
  P_compare_temp <- clean_P(compare.predictors)
  G_compare_temp <- clean_P(compare.group)
  
  if(!is_empty(P_temp_final[P_temp_final %in% P_compare_temp]))
    stop("\n", "Predictor(s) to compare are already in final model, choose new predictors to compare", "\n")
  if(!is_empty(P_temp_final[P_temp_final %in% G_compare_temp]))
    stop("\n", "Predictor(s) to compare are already in final model, choose new predictors to compare", "\n")
  
  P_compare <-
    compare.predictors
  P_group <-
    compare.group
  if(!is_null(compare.group)){
    P_compare <- 1
    P_group <-
      compare.group
  }
  
  nimp <-
    pobj$nimp
  
  Y <-
    c(paste(pobj$Outcome, paste("~")))
  
  nri_res <- fm_step <- auc_stats <- auc_pool_diff <- list()
  D3_stats <- matrix(0, length(P_compare), 3)
  stats_compare <- matrix(0, length(P_compare), 5)
  auc_res <- matrix(0, nimp, 6)
  
  if(boot_auc){
    auc_res <- matrix(0, nimp, 7)
    imp1 <-
      pobj$data[pobj$data[pobj$impvar] == 1, ]
    boot_data <-
      bootstraps(imp1, times = nboot, strata = pobj$Outcome)
  }
  
  n <-
    nrow(pobj$data[pobj$data[pobj$impvar] == 1, ])
  
  Pred_X <-
    array(NA, c(n, 3, nimp))
  
  # Loop k, to pool models in multiply imputed datasets
  for (j in 1:length(P_compare)) {
    
    if(!is_null(P_group)){
      if(any(grepl("[:]", P_group)))
        P_group <-
          gsub(":", "*", P_group)
    }
    if(!is_null(P_compare)){
      if(grepl("[:]", P_compare[j]))
        P_compare[j] <-
          gsub(":", "*", P_compare[j])
    }
    # set regression formula fm
    if(!length(pobj$predictors_final)==0) {
      fm1 <-
        terms(as.formula(paste0(Y, paste0(c(pobj$predictors_final, P_compare[j]), collapse = "+"))))
      fm0 <-
        terms(as.formula(paste0(Y, paste0(c(pobj$predictors_final), collapse = "+"))))
      if(!is_null(P_group))
        fm1 <-
        terms(as.formula(paste0(Y, paste0(c(pobj$predictors_final, P_group), collapse = "+"))))
    } else {
      fm1 <-
        terms(as.formula(paste0(Y, paste0(P_compare[j]), collapse = "+")))
      fm0 <-
        terms(as.formula(paste0(Y, paste0("1"), collapse = "+")))
      if(!is_null(P_group))
        fm1 <-
        terms(as.formula(paste0(Y, paste0(P_group, collapse = "+"))))
    }
    fit1 <- fit0 <- list()
    
    for (i in 1:nimp) {
      data_imp <-
        pobj$data[pobj$data[pobj$impvar] == i, ]
      
      fit1[[i]] <-
        glm(fm1, data = data_imp, x=TRUE, y=TRUE, family = binomial)
      fit0[[i]] <-
        glm(fm0, data = data_imp, x=TRUE, family = binomial)
      Y_imp <-
        data_imp[, pobj$Outcome]
      
      Pred0 <-
        predict(fit0[[i]], type = "response")
      Pred1 <-
        predict(fit1[[i]], type = "response")
      
      # Determine AUC
      roc0 <-
        pROC::roc(Y_imp, Pred0, quiet=TRUE)
      roc1 <-
        pROC::roc(Y_imp, Pred1, quiet=TRUE)
      auc_cov <-
        pROC::cov(roc0, roc1)
      
      auc0 <-
        pROC::roc(Y_imp, Pred0, quiet=TRUE)$auc
      auc0_se <-
        sqrt(var(auc0))
      auc1 <-
        pROC::roc(Y_imp, Pred1, quiet=TRUE)$auc
      auc1_se <-
        sqrt(var(auc1))
      
      auc_se_delong <-
        sqrt(var(auc0) + var(auc1) - 2*auc_cov) # Method DeLong
      
      auc_se_boot <- NA
      if(boot_auc){
        
        # Determine SE AUC with bootstrapping
        boot_diff_auc <- lapply(boot_data$splits, function(x) {
          boot_id <-
            x[[2]]
          boot_data <-
            data_imp[boot_id, ]
          Pred0 <-
            predict(fit0[[i]], type = "response", newdata = boot_data)
          Pred1 <-
            predict(fit1[[i]], type = "response", newdata = boot_data)
          auc0 <-
            pROC::roc(fit0[[i]]$y, Pred0, quiet=TRUE)$auc
          auc1 <-
            pROC::roc(fit1[[i]]$y, Pred1, quiet=TRUE)$auc
          d_auc <-
            auc1-auc0
          return(d_auc)
        })
        auc_se_boot <-
          sd(unlist(boot_diff_auc))
      }
      
      diff_auc <- auc1-auc0
      if(boot_auc)
        auc_res[i, ] <-
        c(auc0, auc0_se,
          auc1, auc1_se,
          diff_auc,
          auc_se_delong, auc_se_boot)
      else
        auc_res[i, ] <-
        c(auc0, auc0_se,
          auc1, auc1_se,
          diff_auc,
          auc_se_delong)
      
      Pred_X[,,i] <-
        cbind(Pred0, Pred1, Y_imp)
    }
    
    # Pool AUC
    auc0_pool <-
      pool_auc(auc_res[, 1], auc_res[, 2], nimp = nimp)
    auc1_pool <-
      pool_auc(auc_res[, 3], auc_res[, 4], nimp = nimp)
    auc_pool <-
      data.frame(rbind(auc0_pool, auc1_pool))
    names(auc_pool) <- c("95% Low", "AUC", "95% Up")
    row.names(auc_pool) <- c("AUC M0", "AUC M1")
    auc_stats[[j]] <-
      auc_pool
    names(auc_stats)[j] <-
      P_compare[j]
    if(!is_null(compare.group))
      names(auc_stats)[j] <- "Group of Predictors"
    
    auc_delong <-
      data.frame(matrix(round(RR_diff_prop(auc_res[, 5], auc_res[, 6]), 5), 1, 4))
    names(auc_delong) <- c("AUC", "SE", "95% CI LO", "95% CI HI")
    row.names(auc_delong) <- "Diff DeLong"
    auc_pool_diff[[j]] <- auc_delong
    if(boot_auc){
      auc_boot <-
        data.frame(matrix(round(RR_diff_prop(auc_res[, 5], auc_res[, 7]), 5), 1, 4))
      names(auc_boot) <- c("AUC", "SE", "95% CI LO", "95% CI HI")
      row.names(auc_boot) <- "Diff Boot"
      auc_pool_diff[[j]] <- rbind(auc_delong, auc_boot)
    }
    names(auc_pool_diff)[j] <-
      P_compare[j]
    if(!is_null(compare.group))
      names(auc_pool_diff)[j] <- "Group of Predictors"
    
    Pred_list <-
      purrr::array_tree(Pred_X, margin = 3)
    
    # Pool NRI
    nri_res[[j]] <-
      pool_reclassification(Pred_list, cutoff=cutoff)
    names(nri_res)[j] <-
      P_compare[j]
    if(!is_null(compare.group))
      names(nri_res)[j] <- "Group of Predictors"
    
    AIC1 <-
      mean(unlist(lapply(fit1, AIC)))
    AIC0 <-
      mean(unlist(lapply(fit0, AIC)))
    
    LogLik1 <-
      mean(-2 * unlist(lapply(fit1, function(x) logLik(x)[1])))
    LogLik0 <-
      mean(-2 * unlist(lapply(fit0, function(x) logLik(x)[1])))
    
    #################### D3 ###################
    
    fit1 <- getfit(fit1)
    m <- length(fit1)
    est1 <- pool(fit1, dfcom = NULL)
    qbar1 <- getqbar(est1)
    
    fit0 <- getfit(fit0)
    
    est0 <- pool(fit0, dfcom = NULL)
    qbar0 <- getqbar(est0)
    k <- length(qbar1) - length(qbar0)
    
    dev1.M <-
      -2 * unlist(lapply(fit1, function(x) logLik(x)[1]))
    dev0.M <-
      -2 * unlist(lapply(fit0, function(x) logLik(x)[1]))
    
    ################# fix.coef ###############
    
    fix.coef_adj <- function(model, beta = NULL)
    {
      est <- tidy(model, effects = "fixed")
      oldcoef <- est$estimate
      names(oldcoef) <- est$term
      #oldcoef <- tidy.coef(model)
      
      beta <- beta[names(oldcoef)]
      data <- model.frame(formula = formula(model), data = model$data)
      mm <- model.matrix(formula(model, fixed.only = TRUE), data = data)
      offset <- as.vector(mm %*% beta)
      uf <- . ~ 1
      upd <- update(model, formula. = uf, data = cbind(data, offset = offset),
                    offset = offset)
      upd
    }
    
    mds1 <- lapply(fit1, fix.coef_adj, beta = qbar1)
    dev1.L <- -2 * unlist(lapply(mds1, function(x) logLik(x)[1]))
    mds0 <- lapply(fit0, fix.coef_adj, beta = qbar0)
    dev0.L <- -2 * unlist(lapply(mds0, function(x) logLik(x)[1]))
    deviances <- list(dev1.M = dev1.M, dev0.M = dev0.M, dev1.L = dev1.L,
                      dev0.L = dev0.L)
    
    dev.M <- mean(dev0.M - dev1.M)
    dev.L <- mean(dev0.L - dev1.L)
    rm <- ((m + 1)/(k * (m - 1))) * (dev.M - dev.L)
    Dm <- dev.L/(k * (1 + rm))
    v <- k * (m - 1)
    if (v > 4) {
      w <- 4 + (v - 4) * ((1 + (1 - 2/v) * (1/rm))^2)
    } else
      w <- v * (1 + 1/k) * ((1 + 1/rm)^2)/2
    pvalue = 1 - pf(Dm, k, w)
    
    D3_res <- c(pvalue, Dm, dev.L)
    D3_stats[j, ] <- D3_res
    AIC_diff <- dev.L - 2*k
    
    stats_compare[j, ] <-
      c(LogLik0, LogLik1, AIC0, AIC1, AIC_diff)
    
    ###############################
    
    # Extract regression formula's
    fm_step[[j]] <-
      paste(Y, paste(attr(fm1, "term.labels"), collapse = " + "))
    names(fm_step)[j] <-
      paste("Test - ", P_compare[j])
    if(!is_null(compare.group))
      names(fm_step)[j] <- paste("Test - ", "Group of Predictors")
  }
  # End j loop
  ##############################################################
  
  D3_stats <-
    round(data.frame(D3_stats), 5)
  row.names(D3_stats) <-
    P_compare
  names(D3_stats) <-
    c("p-values", "D3 statistic", "LRT fixed")
  stats_compare <-
    data.frame(stats_compare)
  row.names(stats_compare) <-
    P_compare
  names(stats_compare) <-
    c("LogLik0", "Loglik1", "AIC0", "AIC1", "AIC_diff")
  if(!is_null(P_group)) {
    row.names(D3_stats) <- row.names(stats_compare) <- NULL
    P_compare <- NULL
  }
  
  fm_null <-
    formula(fm0)
  if(length(pobj$predictors_final)==0)
    fm_null <-
    as.formula(paste0(pobj$Outcome, paste0("~ 1")))
  
  compobj <-
    list(D3_stats=D3_stats, stats_compare=stats_compare,
         NRI=nri_res, cutoff=cutoff, AUC = auc_stats, 
         AUC_diff = auc_pool_diff, formula_test=fm_step, 
         null_model = fm_null, compare_predictors = P_compare, 
         compare_group = compare.group)
  compobj
}