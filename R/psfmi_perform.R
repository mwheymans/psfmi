#' Evaluate performance of logistic regression models in Multiply Imputed datasets
#'
#' \code{psfmi_perform} Evaluate Performance of logistic regression models selected with
#'  the \code{psfmi_lr} function of the \code{psfmi} package.
#'
#' @param pobj An object of class \code{smodsmi} (selected models in multiply imputed datasets), 
#'  produced by a previous call to \code{psfmi_lr}.
#' @param data_orig dataframe of original dataset that contains missing data for method boot_MI 
#' @param nboot The number of bootstrap resamples, default is 10.
#' @param int_val If TRUE internal validation is conducted in multiply imputed datasets. 
#'  See \code{method} for methods that can be used.
#' @param method Methods for internal validation in multiply imputed datasets.
#'  Choose MI_boot for bootstrapping in each imputed dataset and boot_MI for multiple 
#'  imputation in each bootstrap sample. To use the second method data_orig has to be specified.
#'  The first method is faster. See details for more information.
#' @param nimp_boot_MI Numerical scalar. Number of imputed datasets for method boot_MI.
#'  Default is 5. When not defined, the number of multiply imputed datasets is used of the  
#'  previous call to the function code{psfmi_lr}.
#' @param p.crit A numerical scalar. P-value selection criterium used for backward selection
#'  during internal validation. When set at 1, pooling and internal validation is done without 
#'  backward selection. 
#' @param mice_method The Multiple Imputation method used for each predictor with missing values. 
#'  For Multiple Imputation the \code{mice} package is used. See that package for more information.      
#' @param mice_niter Numerical scalar. Default is 10. The number of iterations in Multiple Imputation. 
#'  See the \code{mice} package for more information.       
#' @param mice_seed Numerical scalar. Default is random number generator initializeb by computer 
#'  via set.seed().      
#' @param predictorMatrix A numeric matrix of nrow(data) rows and ncol(data) columns, containing 0/1 
#'  data specifying the imputation models used to impute the predictors with missing data. Default 
#'  is that each variable is used to impute other variables. See the \code{mice} package 
#'  for more information.     
#' @param cal.plot If TRUE a calibration plot is generated. Default is FALSE. Can be used in combination
#'  with int_val = FALSE.
#' @param plot.indiv If TRUE calibration plots for each separate imputed dataset are generated, 
#'  otherwise all calibration plots are plotted in one figure.       
#' @param groups_cal A numerical scalar. Number of groups used on the calibration plot. 
#'  Default is 10. If the range of predicted probabilities is too low 5 groups can be
#'  chosen.       
#'  
#' @details For internal validation two methods can be used, MI_boot and boot_MI. MI_boot draws
#'  for each bootstrap step the same cases in all imputed datasets. With boot_MI first bootstrap samples
#'  are drawn from the original dataset with missing values and than multiple imputation is
#'  applied. For multiple imputation the \code{mice} function from the \code{mice} package is used. 
#'  It is recommended to use a minumum of 100 bootstrap samples, which may take some time. The method
#'  boot_MI is more time consuming than MI_boot.      
#'                 
#'@return A \code{psfmi_perform} object from which the following objects can be extracted: \code{res_boot}, 
#'  result of pooled performance (in multiply imputed datasets) at each bootstrap step of ROC app (pooled 
#'  ROC), ROC test (pooled ROC after bootstrap model is applied in original multiply imputed datasets), 
#'  same for R2 app (Nagelkerke's R2),  R2 test, Brier app and Brier test. Information is also provided 
#'  about testing the Calibration slope at each bootstrap step as interc test and 
#'  Slope test. The performance measures are pooled by a call to the function \code{pool_performance}. Another
#'  object that can be extracted is \code{intval}, with information of the AUC, R2, Brier score and 
#'  Calibration slope averaged over the bootstrap samples, in terms of: Orig (original datasets), 
#'  Apparent (models applied in bootstrap samples), Test (bootstrap models are applied in original datasets),
#'  Optimism (difference between apparent and test) and Corrected (original corrected for optimism). 
#' 
#' @references Heymans MW, van Buuren S, Knol DL, van Mechelen W, de Vet HC. Variable selection under 
#'  multiple imputation using the bootstrap in a prognostic study. BMC Med Res Methodol. 2007(13);7:33.
#' @references F. Harrell. Regression Modeling Strategies. With Applications to
#'  Linear Models, Logistic and Ordinal Regression, and Survival Analysis (2nd edition). Springer,
#'  New York, NY, 2015.
#' @references Van Buuren S. (2018). Flexible Imputation of Missing Data. 2nd Edition. Chapman &
#'  Hall/CRC Interdisciplinary Statistics. Boca Raton.
#' @references Harel, O. (2009). The estimation of R2 and adjusted R2 in
#'  incomplete data sets using multiple imputation. Journal of Applied Statistics,
#'  36(10), 1109-1118.
#' @references Musoro JZ, Zwinderman AH, Puhan MA, ter Riet G, Geskus RB. Validation of prediction 
#'  models based on lasso regression with multiply imputed data. BMC Med Res Methodol. 2014;14:116. 
#' @references Wahl S, Boulesteix AL, Zierer A, Thorand B, van de Wiel MA. Assessment of 
#'  predictive performance in incomplete data by combining internal validation and multiple 
#'  imputation. BMC Med Res Methodol. 2016;16(1):144. 
#' @references EW. Steyerberg (2019). Clinical Prediction MOdels. A Practical Approach 
#'  to Development, Validation, and Updating (2nd edition). Springer Nature Switzerland AG.
#'  
#' @references http://missingdatasolutions.rbind.io/
#' 
#' @examples
#' res_psfmi <- psfmi_lr(data=lbpmilr, nimp=5, impvar="Impnr", Outcome="Chronic", 
#'   predictors=c("Gender", "Pain","Tampascale","Smoking","Function", "Radiation", 
#'   "Age"), p.crit = 1, method="D1")
#'  
#' res_val <- psfmi_perform(res_psfmi, int_val = FALSE, p.crit=1, cal.plot=TRUE, 
#' plot.indiv=FALSE)
#' res_val  
#' 
#' res_psfmi <- psfmi_lr(data=lbpmilr, nimp=5, impvar="Impnr", Outcome="Chronic", 
#'   predictors=c("Gender", "Pain","Tampascale","Smoking","Function", "Radiation", 
#'   "Age"), cat.predictors = "Carrying", keep.predictors = "Function", 
#'   p.crit = 0.157, method="D1")   
#'  
#' res_val <- psfmi_perform(res_psfmi, int_val = TRUE, method = "MI_boot", 
#'   nboot = 10, p.crit=0.157,
#'   cal.plot=FALSE, plot.indiv=FALSE)
#'   res_val$res_boot
#'   res_val$intval
#'  
#' @export
psfmi_perform <- function(pobj, data_orig = NULL, nboot = 10, int_val = FALSE, method = NULL, 
                      nimp_boot_MI = 5, p.crit = 1, mice_method = NULL, mice_niter = 10, 
                      mice_seed = NA, predictorMatrix=NULL, cal.plot=FALSE, plot.indiv=FALSE, 
                      groups_cal=10)
{
##############################
  call <- match.call()
if(class(pobj)!="smodsmi")
  stop("\n", "Object should be of type smodsmi", "\n")
if(pobj$method=="D2")
  stop("\n", "Choose method D1, D3 or MPR for pooling or variable selection first,", "\n",
         "Method D2 can become unstable during bootstrap validation", "\n")
if(pobj$model_type=="survival")
  stop("\n", "Function only available for models of type binomial", "\n")
if(!is.null(pobj$random.eff))
    stop("\n", "Function only available for regression models without random effects", "\n")
  
if(int_val==TRUE) {
  cal.plot==FALSE
  message("\n", "Calibration plot not made, only when int.val = FALSE")
}
############################################################################### Method MI_boot

if(int_val==TRUE){

  if(is.null(nboot))
    stop("\n", "Number of bootstrap samples not defined")
  if(nboot<1)
    stop("\n", "Increase number of bootstrap samples")
  if(is.null(method))
    stop("\n", "Validation method not defined, choose boot_MI or MI_boot")

if(method=="MI_boot") {

  if(is.null(nimp_boot_MI))
    nimp_boot_MI <- pobj$nimp
  if (nimp_boot_MI < 2)
    stop("\n", "Number of imputed datasets for MI_boot must be > 1", "\n")

imp1 <- pobj$data[pobj$data[pobj$impvar] == 1, ]
boot_data <- bootstraps(imp1, times = nboot)

boot_seq <- as.list(1:nboot)

opt_boot <- mapply(function(x, y) {
  message("\n", "Boot ", y)
  boot_id <- x[[2]]

  # Split original imputed data and select
  # same bootstrap cases in imputed datasets
  split_orig_imp <- pobj$data %>%
    group_by(pobj$data[,pobj$impvar])
  split_app_data <- group_split(split_orig_imp)
  # select in each original imputed datasets the same bootstrap id's
  boot_app_data <- lapply(split_app_data, function(x) {
    x <- as.data.frame(x)
    app_data_select <- x[boot_id, ]
  })
  app_data <- do.call("rbind", boot_app_data)

  pool_model_lr <- psfmi_lr(data=app_data, nimp=pobj$nimp, impvar = pobj$impvar,
                            Outcome = pobj$Outcome, predictors = pobj$predictors,
                            p.crit = p.crit, cat.predictors = pobj$cat.predictors,
                            spline.predictors = pobj$spline.predictors,
                            int.predictors = pobj$int.predictors, keep.predictors = pobj$keep.predictors,
                            knots = pobj$knots, method = pobj$method, print.method = pobj$print.method)

  predictors_selected <- pool_model_lr$predictors_in[nrow(pool_model_lr$predictors_in), ]
  names_predictors_selected <- pool_model_lr$predictors_initial

  Y <- c(paste(pobj$Outcome, paste("~")))
  if(is_empty(pool_model_lr$predictors_final)) {
    #message("Model is empty")
    pool_model_lr$predictors_final <- 1
    fm <- as.formula(paste(Y, paste(pool_model_lr$predictors_final, collapse = "+")))
    lp_app_pooled <- 1
  } else {
    fm <- as.formula(paste(Y, paste(pool_model_lr$predictors_final, collapse = "+")))
    lp_app_pooled <- c(pool_model_lr$RR_Model[[length(pool_model_lr$RR_Model)]][, 1])
  }

  # Obtain apparent pooled performance measures
  perform_app <- pool_performance(data=app_data, nimp = pobj$nimp,
                                  impvar=pobj$impvar, Outcome = pobj$Outcome, 
                                  predictors = pool_model_lr$predictors_final, cal.plot = FALSE)

  # Test apparent LP in (imputed) original data
  coef_slope_test <- list()
  coef.f <- se.f <- pred.i <- rsq.i <- brier.i <- list()
  roc.f.i <- se.roc.f.i <- se.roc.f.i.logit <- list()

  for(i in 1:pobj$nimp){
    data_test <- pobj$data[pobj$data[pobj$impvar] == i, ]

    fit <- glm(fm, data = data_test, family = binomial)

    lp_test <- model.matrix(fit) %*% lp_app_pooled
    fit_test <- glm(fit$y ~ lp_test, family = binomial)
    coef_fit_test <- coef(fit_test)

    if(is_empty(pool_model_lr$predictors_final)) {
      coef_fit_test <- replace_na(coef_fit_test, 1)
    }
    coef_slope_test[[i]] <- coef_fit_test

    test_prob <- c(1/(1+exp(-lp_test)))

    # ROC/AUC
    roc.f.i[[i]] <- roc(fit$y, test_prob, quiet = TRUE)$auc
    se.roc.f.i[[i]] <- sqrt(pROC::var(roc.f.i[[i]]))
    se.roc.f.i.logit[[i]] <- sqrt(pROC::var(roc.f.i[[i]])) /
      (roc.f.i[[i]]*(1-roc.f.i[[i]]))

    # Nagelkerke R-squared
    f.full <- -2 * logLik(fit_test)
    f.base <- -2 * logLik(update(fit_test, ~ 1))
    n <- fit_test$df.null
    rsq.nagel <- (1 - exp((f.full - f.base)/n))/
      (1 - exp(-f.base/n))
    rsq.i[[i]] <- rsq.nagel

    # Brier and scaled Brier score
    brier <- function(obs, pred) {
      mean((obs - pred)^2)
    }
    brier.i[[i]] <- brier(fit$y, test_prob)
  }

  # ROC/AUC
  # Rubin's Rules on logit transformation ROC curve and SE
  est.roc.logit <- log(unlist(roc.f.i)/
                         (1-unlist(roc.f.i)))
  se.roc.logit <- unlist(se.roc.f.i.logit)

  # Pooling
  p.roc.logit <- mean(est.roc.logit)
  # within variance
  p.se.roc.logit <- mean(se.roc.logit)
  # between variance
  b.roc.logit <- var(est.roc.logit)
  # total variance
  tv.roc.logit <- p.se.roc.logit +
    ((1 + (1/pobj$nimp)) * b.roc.logit)
  se.t.roc.logit <- sqrt(tv.roc.logit)

  # Backtransform
  inv.roc <- exp(p.roc.logit) /
    (1 + exp(p.roc.logit))
  inv.roc.u <- exp(p.roc.logit + (1.96*se.t.roc.logit)) /
    (1 + exp(p.roc.logit + (1.96*se.t.roc.logit)))
  inv.roc.l <- exp(p.roc.logit - (1.96*se.t.roc.logit)) /
    (1 + exp(p.roc.logit - (1.96*se.t.roc.logit)))

  roc_test <- matrix(c(inv.roc.l, inv.roc, inv.roc.u),
                     1, 3, byrow = T)
  dimnames(roc_test) <- list(c("ROC (logit)"),
                             c("95% Low", "ROC", "95% Up"))

  # Pooling R square
  # Fisher z Transformation
  z.rsq <- atanh(unlist(rsq.i))
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

  brier_pool_test <- mean(unlist(brier.i))

  lp_test_pooled <- colMeans(do.call("rbind", coef_slope_test), na.rm = TRUE)
  # End pooling performance measures in multiply imputed data

  ROC_app <- perform_app$ROC_pooled[2]
  ROC_test <- roc_test[2]
  Slope_test <- lp_test_pooled
  R2_app <- perform_app$R2_pooled
  Brier_app <- perform_app$Brier_pooled
  Brier_test <- brier_pool_test

  opt_perform <- list(c(ROC_app, ROC_test, R2_app, R2_test,
                   Brier_app, Brier_test, Slope_test),
                   predictors_selected, names_predictors_selected)
  return(opt_perform)

}, x = boot_data$splits, y=boot_seq, SIMPLIFY = FALSE)

predictors_selected <- data.frame(do.call("rbind", lapply(opt_boot, function(x) x[[2]])))
colnames(predictors_selected) <- opt_boot[[1]][[3]]
row.names(predictors_selected) <- paste("Boot", 1:nboot)

res_boot <- data.frame(do.call("rbind", lapply(opt_boot, function(x) x[[1]])))
colnames(res_boot) <- c("ROC app", "ROC test", "R2 app", "R2 test",
                        "Brier app", "Brier test", "interc test", "Slope test")
row.names(res_boot) <- paste("Boot", 1:nboot)

roc_opt <- res_boot[,1] - res_boot[,2]
r2_opt <- res_boot[,3] - res_boot[,4]
brier_opt <- res_boot[,5] - res_boot[,6]

res_boot_m <- colMeans(data.frame(res_boot, roc_opt, r2_opt, brier_opt), na.rm=TRUE)

# Perform original model in multiply imputed original data
Y_orig <- c(paste(pobj$Outcome, paste("~")))
if(is_empty(pobj$predictors_final)) {
  pobj$predictors_final <- 1
  fm_orig <- as.formula(paste(Y_orig, paste(pobj$predictors_final, collapse = "+")))
} else {
fm_orig <- as.formula(paste(Y_orig, paste(pobj$predictors_final, collapse = "+")))
}
perform_mi_orig <- pool_performance(data=pobj$data, nimp = pobj$nimp,
                                    impvar=pobj$impvar, Outcome = pobj$Outcome,
                                    predictors = pobj$predictors_final, cal.plot = FALSE)

ROC_orig <- perform_mi_orig$ROC_pooled[2]
R2_orig <-perform_mi_orig$R2_pooled
Brier_orig <-perform_mi_orig$Brier_pooled

ROC_corr <- ROC_orig - res_boot_m[9]
ROC_val <- c(ROC_orig, res_boot_m[1], res_boot_m[2], res_boot_m[9], ROC_corr)

R2_corr <- R2_orig - res_boot_m[10]
R2_val <- c(R2_orig, res_boot_m[3], res_boot_m[4], res_boot_m[10], R2_corr)

Brier_corr <- Brier_orig - res_boot_m[11]
Brier_val <- c(Brier_orig, res_boot_m[5], res_boot_m[6], res_boot_m[11], Brier_corr)

Slope_corr <- res_boot_m[8]
Slope_val <- c(1, 1, res_boot_m[8], 1-res_boot_m[8], res_boot_m[8])

val_res <- as.data.frame(matrix(c(ROC_val, R2_val, Brier_val, Slope_val), 4, 5, byrow = T))

#print(predictors_selected)

colnames(val_res) <- c("Orig", "Apparent", "Test", "Optimism", "Corrected")
row.names(val_res) <- c("AUC", "R2", "Brier", "Slope")

}
################################################################################ Method boot_MI
if (method=="boot_MI"){

if(is.null(data_orig))
    stop("Include original data with missing values for the boot_MI method")
if(is.null(nimp_boot_MI))
  nimp_boot_MI <- pobj$nimp

data_orig <- data_orig
if(is.null(predictorMatrix))
  predictorMatrix <- mice::make.predictorMatrix(data_orig)

# bootstrap from original data with missings included
boot_data_orig <- bootstraps(data_orig, times = nboot)

boot_seq <- as.list(1:nboot)
# use bootstrap data as input
# Start bootstrap run
opt_boot <- mapply(function(x, y) {
  message("\n", "Boot ", y)

  # x is bootstrap sample with missing data
  x <- as.data.frame(x)

  # Multiply Impute bootstrap sample with mice
  boot_data_imp <- mice(x, m=nimp_boot_MI, method = mice_method,
          maxit = mice_niter, predictorMatrix = predictorMatrix, seed = mice_seed, printFlag = FALSE)
  boot_data_compl <- complete(boot_data_imp, action = "long", include = FALSE)

  # Pool model in completed datasets (set p.crit = 1 to exclude bws)
  pool_model_lr <- psfmi_lr(data=boot_data_compl, nimp=nimp_boot_MI, impvar = ".imp",
      Outcome = pobj$Outcome, predictors = pobj$predictors,
      p.crit = p.crit, cat.predictors = pobj$cat.predictors,
      spline.predictors = pobj$spline.predictors,
      int.predictors = pobj$int.predictors, keep.predictors = pobj$keep.predictors,
      knots = pobj$knots, method = pobj$method, print.method = pobj$print.method)

  predictors_selected <- pool_model_lr$predictors_in[nrow(pool_model_lr$predictors_in), ]
  names_predictors_selected <- pool_model_lr$predictors_initial

  # Extract apparent pooled LP (LP in each bootstrap sample)
  Y <- c(paste(pobj$Outcome, paste("~")))
  if(is_empty(pool_model_lr$predictors_final)) {
    #message("Model is empty")
    pool_model_lr$predictors_final <- 1
    fm <- as.formula(paste(Y, paste(pool_model_lr$predictors_final, collapse = "+")))
    lp_app_pooled <- 1
  } else {
    fm <- as.formula(paste(Y, paste(pool_model_lr$predictors_final, collapse = "+")))
    lp_app_pooled <- c(pool_model_lr$RR_Model[[length(pool_model_lr$RR_Model)]][, 1])
  }

  # Obtain apparent pooled performance measures
  perform_app <- pool_performance(data=boot_data_compl, nimp = nimp_boot_MI,
            impvar=".imp", Outcome = pobj$Outcome,
            predictors = pool_model_lr$predictors_final, cal.plot = FALSE)

  # Test apparent LP in original Multiply Imputed data
  coef_slope_test <- list()
  coef.f <- se.f <- pred.i <- rsq.i <- brier.i <- list()
  roc.f.i <- se.roc.f.i <- se.roc.f.i.logit <- list()

  for(i in 1:pobj$nimp){
    data_test <- pobj$data[pobj$data[pobj$impvar] == i, ]
    fit <- glm(fm, data = data_test, family = binomial)

    lp_test <- model.matrix(fit) %*% lp_app_pooled
    fit_test <- glm(fit$y ~ lp_test, family = binomial)
    coef_fit_test <- coef(fit_test)

    if(is_empty(pool_model_lr$predictors_final)) {
      coef_fit_test <- replace_na(coef_fit_test, 1)
    }
    coef_slope_test[[i]] <- coef_fit_test

    test_prob <- c(1/(1+exp(-lp_test)))

    # ROC/AUC
    roc.f.i[[i]] <- roc(fit$y, test_prob, quiet = TRUE)$auc
    se.roc.f.i[[i]] <- sqrt(pROC::var(roc.f.i[[i]]))
    se.roc.f.i.logit[[i]] <- sqrt(pROC::var(roc.f.i[[i]])) /
      (roc.f.i[[i]]*(1-roc.f.i[[i]]))

    # Nagelkerke R-squared
    f.full <- -2 * logLik(fit_test)
    f.base <- -2 * logLik(update(fit_test, ~ 1))
    n <- fit_test$df.null
    rsq.nagel <- (1 - exp((f.full - f.base)/n))/
      (1 - exp(-f.base/n))
    rsq.i[[i]] <- rsq.nagel

    # Brier and scaled Brier score
    brier <- function(obs, pred) {
      mean((obs - pred)^2)
    }
    brier.i[[i]] <- brier(fit$y, test_prob)

  }
  # End pooling performance measures in multiply imputed data

  # ROC/AUC
  # Rubin's Rules on logit transformation ROC curve and SE
  est.roc.logit <- log(unlist(roc.f.i)/
                         (1-unlist(roc.f.i)))
  se.roc.logit <- unlist(se.roc.f.i.logit)

  # Pooling
  p.roc.logit <- mean(est.roc.logit)
  # within variance
  p.se.roc.logit <- mean(se.roc.logit)
  # between variance
  b.roc.logit <- var(est.roc.logit)
  # total variance
  tv.roc.logit <- p.se.roc.logit +
    ((1 + (1/pobj$nimp)) * b.roc.logit)
  se.t.roc.logit <- sqrt(tv.roc.logit)

  # Backtransform
  inv.roc <- exp(p.roc.logit) /
    (1 + exp(p.roc.logit))
  inv.roc.u <- exp(p.roc.logit + (1.96*se.t.roc.logit)) /
    (1 + exp(p.roc.logit + (1.96*se.t.roc.logit)))
  inv.roc.l <- exp(p.roc.logit - (1.96*se.t.roc.logit)) /
    (1 + exp(p.roc.logit - (1.96*se.t.roc.logit)))

  roc_test <- matrix(c(inv.roc.l, inv.roc, inv.roc.u),
                      1, 3, byrow = T)
  dimnames(roc_test) <- list(c("ROC (logit)"),
                              c("95% Low", "ROC", "95% Up"))

  # Pooling R square
  # Fisher z Transformation
  z.rsq <- atanh(unlist(rsq.i))
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

  brier_pool_test <- mean(unlist(brier.i))

  lp_test_pooled <- colMeans(do.call("rbind", coef_slope_test), na.rm = TRUE)

  # End pooled test performance measures

  ROC_app <- perform_app$ROC_pooled[2]
  ROC_test <- roc_test[2]
  Slope_test <- lp_test_pooled
  R2_app <- perform_app$R2_pooled
  Brier_app <- perform_app$Brier_pooled
  Brier_test <- brier_pool_test

  opt_perform <- list(c(ROC_app, ROC_test, R2_app, R2_test,
                        Brier_app, Brier_test, Slope_test),
                      predictors_selected, names_predictors_selected)
  return(opt_perform)

  }, x = boot_data_orig$splits, y=boot_seq, SIMPLIFY = FALSE)

predictors_selected <- data.frame(do.call("rbind", lapply(opt_boot, function(x) x[[2]])))
colnames(predictors_selected) <- opt_boot[[1]][[3]]
row.names(predictors_selected) <- paste("Boot", 1:nboot)

res_boot <- data.frame(do.call("rbind", lapply(opt_boot, function(x) x[[1]])))
colnames(res_boot) <- c("ROC app", "ROC test", "R2 app", "R2 test",
                        "Brier app", "Brier test", "interc test", "Slope test")
row.names(res_boot) <- paste("Boot", 1:nboot)

roc_opt <- res_boot[,1] - res_boot[,2]
r2_opt <- res_boot[,3] - res_boot[,4]
brier_opt <- res_boot[,5] - res_boot[,6]

res_boot_m <- colMeans(data.frame(res_boot, roc_opt, r2_opt, brier_opt), na.rm=TRUE)

# Perform original model in multiply imputed original data
Y_orig <- c(paste(pobj$Outcome, paste("~")))
if(is_empty(pobj$predictors_final)) {
  pobj$predictors_final <- 1
  fm_orig <- as.formula(paste(Y_orig, paste(pobj$predictors_final, collapse = "+")))
} else {
  fm_orig <- as.formula(paste(Y_orig, paste(pobj$predictors_final, collapse = "+")))
}
perform_mi_orig <- pool_performance(data=pobj$data, nimp = pobj$nimp,
                                    impvar=pobj$impvar, Outcome = pobj$Outcome,
                                    predictors = pobj$predictors_final, cal.plot = FALSE)

ROC_orig <- perform_mi_orig$ROC_pooled[2]
R2_orig <-perform_mi_orig$R2_pooled
Brier_orig <-perform_mi_orig$Brier_pooled

ROC_corr <- ROC_orig - res_boot_m[9]
ROC_val <- c(ROC_orig, res_boot_m[1], res_boot_m[2], res_boot_m[9], ROC_corr)

R2_corr <- R2_orig - res_boot_m[10]
R2_val <- c(R2_orig, res_boot_m[3], res_boot_m[4], res_boot_m[10], R2_corr)

Brier_corr <- Brier_orig - res_boot_m[11]
Brier_val <- c(Brier_orig, res_boot_m[5], res_boot_m[6], res_boot_m[11], Brier_corr)

Slope_corr <- res_boot_m[8]
Slope_val <- c(1, 1, res_boot_m[8], 1-res_boot_m[8], res_boot_m[8])

val_res <- as.data.frame(matrix(c(ROC_val, R2_val, Brier_val, Slope_val), 4, 5, byrow = T))

colnames(val_res) <- c("Orig", "Apparent", "Test", "Optimism", "Corrected")
row.names(val_res) <- c("AUC", "R2", "Brier", "Slope")
}

pobjval <- list(res_boot = res_boot, intval = val_res, 
                predictors_in=predictors_selected, nboot = nboot)
return(pobjval)
}
if(!int_val){
  message("\n", "Performance measures pooled in original multiply imputed datasets", "\n")
  Y <- c(paste(pobj$Outcome, paste("~")))
  if(is_empty(pobj$predictors_final)) {
    pobj$predictors_final <- 1
    fm <- as.formula(paste(Y, paste(pobj$predictors_final, collapse = "+")))
  } else {
    fm <- as.formula(paste(Y, paste(pobj$predictors_final, collapse = "+")))
  }

  perform_mi_orig <- pool_performance(data=pobj$data, nimp = pobj$nimp,
                                      impvar=pobj$impvar, Outcome = pobj$Outcome,
                                      predictors = pobj$predictors_final, cal.plot=cal.plot, 
                                      plot.indiv=plot.indiv,
                                      groups_cal = groups_cal)
}
return(perform_mi_orig)
}
