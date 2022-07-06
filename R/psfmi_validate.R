#' Internal validation and performance of logistic prediction models across Multiply Imputed datasets
#'
#' \code{psfmi_validate} Evaluate Performance of logistic regression models selected with
#'  the \code{psfmi_lr} function of the \code{psfmi} package by using cross-validation
#'  or bootstrapping.
#'
#' @param pobj An object of class \code{pmods} (pooled models), produced by a previous
#'  call to \code{psfmi_lr}.
#' @param val_method Method for internal validation. MI_boot for first Multiple Imputation and
#'  than bootstrapping in each imputed dataset and boot_MI for first bootstrapping and than
#'  multiple imputation in each bootstrap sample, and cv_MI, cv_MI_RR and
#'  MI_cv_naive for the combinations of cross-validation and multiple imputation.
#'  To use cv_MI, cv_MI_RR and boot_MI, data_orig has to be specified. See details for more information.
#' @param data_orig dataframe of original dataset that contains missing data for methods
#'  cv_MI, cv_MI_RR and boot_MI.
#' @param int_val If TRUE internal validation is conducted using bootstrapping or cross-validation.
#'  Default is TRUE. If FALSE only apparent performance measures are calculated.
#' @param nboot The number of bootstrap resamples, default is 10. Used for methods boot_MI and MI_boot.
#' @param folds The number of folds, default is 3. Used for methods cv_MI, cv_MI_RR and MI_cv_naive.
#' @param nimp_cv Numerical scalar. Number of (multiple) imputation runs for method cv_MI.
#' @param nimp_mice Numerical scalar. Number of imputed datasets for method cv_MI_RR and boot_MI.
#'  When not defined, the number of multiply imputed datasets is used of the
#'  previous call to the function \code{psfmi_lr}.
#' @param p.crit A numerical scalar. P-value selection criterium used for backward or forward
#'  selection during validation. When set at 1, pooling and internal validation is done without
#'  backward selection.
#' @param BW Only used for methods cv_MI, cv_MI_RR and MI_cv_naive. If TRUE backward selection is
#'  conducted within cross-validation. Default is FALSE.
#' @param direction Can be used together with val_methods boot_MI and MI_boot. The direction of
#'  predictor selection, "BW" is for backward selection and "FW" for forward selection.
#' @param cv_naive_appt Can be used in combination with val_method MI_cv_naive. Default is TRUE for
#'  showing the cross-validation apparent (train) and test results. Set to FALSE to only give test results.
#' @param cal.plot If TRUE a calibration plot is generated. Default is FALSE. Can be used in combination
#'  with int_val = FALSE.
#' @param plot.method If "mean" one calibration plot is generated, first taking the 
#'   mean of the linear predictor across the multiply imputed datasets (default), if 
#'   "individual" the calibration plot of each imputed dataset is plotted, 
#'   if "overlay" calibration plots from each imputed datasets are plotted in one figure.       
#' @param groups_cal A numerical scalar. Number of groups used on the calibration plot and. 
#'  for the Hosmer and Lemeshow test. Default is 10. If the range of predicted probabilities. 
#'  is low, less than 10 groups can be chosen, but not < 3. 
#' @param miceImp Wrapper function around the \code{mice} function.
#' @param ...  Arguments as predictorMatrix, seed, maxit, etc that can be adjusted for
#'  the \code{mice} function. To be used in combination with validation methods cv_MI,
#'  cv_MI_RR and MI_boot. For method cv_MI the number of imputed datasets is fixed at 1 and cannot
#'  be changed.
#'
#' @details For internal validation five methods can be used, cv_MI, cv_MI_RR, MI_cv_naive,
#'  MI_boot and boot_MI. Method cv_MI uses imputation within each cross-validation fold definition.
#'  By repeating this in several imputation runs, multiply imputed datasets are generated. Method
#'  cv_MI_RR uses multiple imputation within the cross-validation definition. MI_cv_naive, applies
#'  cross-validation within each imputed dataset. MI_boot draws for each bootstrap step the same
#'  cases in all imputed datasets. With boot_MI first bootstrap samples are drawn from the original
#'  dataset with missing values and than multiple imputation is applied. For multiple imputation
#'  the \code{mice} function from the \code{mice} package is used. It is recommended to use a minumum
#'  of 100 imputation runs for method cv_MI or 100 bootstrap samples for method boot_MI or MI_boot.
#'  Methods cv_MI, cv_MI_RR and MI_cv_naive can be combined with backward selection during
#'  cross-validation and with methods boot_MI and MI_boot, backward and forward selection can
#'  be used. For methods cv_MI and cv_MI_RR the outcome in the original dataset has to be complete.
#'
#'@return A \code{psfmi_perform} object from which the following objects can be extracted: \code{res_boot},
#'  result of pooled performance (in multiply imputed datasets) at each bootstrap step of ROC app (pooled
#'  ROC), ROC test (pooled ROC after bootstrap model is applied in original multiply imputed datasets),
#'  same for R2 app (Nagelkerke's R2), R2 test, Scaled Brier app and Scaled Brier test. Information is also provided
#'  about testing the Calibration slope at each bootstrap step as interc test and Slope test.
#'  The performance measures are pooled by a call to the function \code{pool_performance}. Another
#'  object that can be extracted is \code{intval}, with information of the AUC, R2, Scaled Brier score and
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
#' @section Vignettes:
#' 
#' \itemize{
#' \item  \href{https://mwheymans.github.io/psfmi/articles/cv_MI.html}{MI and Cross-validation - Method cv_MI}
#' \item  \href{https://mwheymans.github.io/psfmi/articles/cv_MI_RR.html}{MI and Cross-validation - Method cv_MI_RR}
#' \item  \href{https://mwheymans.github.io/psfmi/articles/MI_cv_naive.html}{MI and Cross-validation - Method MI_cv_naive} 
#' \item  \href{https://mwheymans.github.io/psfmi/articles/boot_MI.html}{MI and Bootstrapping - Method boot_MI}
#' \item  \href{https://mwheymans.github.io/psfmi/articles/MI_boot.html}{MI and Bootstrapping - Method MI_boot}
#' }
#' 
#' @author Martijn Heymans, 2020
#'
#' @examples
#' pool_lr <- psfmi_lr(data=lbpmilr, formula = Chronic ~ Pain + JobDemands + rcs(Tampascale, 3) +
#'            factor(Satisfaction) + Smoking, p.crit = 1, direction="FW",
#'            nimp=5, impvar="Impnr", method="D1")
#'            
#' pool_lr$RR_model
#'
#' res_perf <- psfmi_validate(pool_lr, val_method = "cv_MI", data_orig = lbp_orig, folds=3,
#'             nimp_cv = 2, p.crit=0.05, BW=TRUE, miceImp = miceImp, printFlag = FALSE)
#'             
#' res_perf
#'
#'\dontrun{
#'  set.seed(200)
#'   res_val <- psfmi_validate(pobj, val_method = "boot_MI", data_orig = lbp_orig, nboot = 5,
#'   p.crit=0.05, BW=TRUE, miceImp = miceImp, nimp_mice = 5, printFlag = FALSE, direction = "FW")
#'   
#'   res_val$stats_val
#'}
#'
#' @export
psfmi_validate <- function(pobj, 
                          val_method = NULL, 
                          data_orig = NULL, 
                          int_val = TRUE, 
                          nboot = 10,
                          folds=3, 
                          nimp_cv = 5, 
                          nimp_mice = 5, 
                          p.crit = 1, 
                          BW = FALSE,
                          direction = NULL, 
                          cv_naive_appt=FALSE,
                          cal.plot=FALSE, 
                          plot.method="mean", 
                          groups_cal=5, 
                          miceImp, 
                          ...)
{
  ##############################
  #General Settings
  call <- match.call()
  
  if(!inherits(pobj, "pmods"))
    stop("\n", "Object should be of type pmods", "\n")
  if(pobj$model_type!="binomial")
    stop("\n", "Methods only available for models of type binomial", "\n")
  if(is_empty(pobj$predictors_final))
    stop("\n", "Model is empty. Cannot validate empty model", "\n")
  if(is.null(val_method))
    stop("\n", "Validation method not defined, choose boot_MI, MI_boot, cv_MI, cv_MI_RR or MI_cv_naive")
  
  if(!int_val){
    message("\n", "No validation - Apparent performance", "\n")
    Y <- c(paste(pobj$Outcome, paste("~")))
    if(is_empty(pobj$predictors_final)) {
      pobj$predictors_final <- 1
      fm <- as.formula(paste(Y, paste(pobj$predictors_final, collapse = "+")))
    } else {
      fm <- as.formula(paste(Y, paste(pobj$predictors_final, collapse = "+")))
    }
    
    perform_mi_orig <- pool_performance(data=pobj$data, nimp = pobj$nimp,
                                        impvar=pobj$impvar,
                                        formula = fm,
                                        cal.plot=cal.plot,
                                        plot.method=plot.method,
                                        groups_cal=groups_cal)
    return(perform_mi_orig)
  }
  if(int_val){
    
    # Specific for cv_MI methods
    if(val_method=="cv_MI" | val_method=="cv_MI_RR") {
      if(BW==FALSE & p.crit!=1)
        stop("\n", "If BW == FALSE, p.crit must be 1", "\n")
      if(BW & p.crit==1)
        stop("\n", "If BW == TRUE, p.crit must be < 1","\n")
      if(is_empty(data_orig))
        stop("\n", "data_orig not defined", "\n")
      if(any(is.na(pobj$Outcome)))
        stop("\n", "Outcome variable contains missing data,
             cv_MI or cv_MI_RR can only be used when outcome is complete", "\n")
      if(val_method=="cv_MI"){
        if(nimp_cv==10){
          message("\n", "Recommended to increase nimp_cv to 100", "\n")
        }
        pobjcv <- cv_MI(pobj, data_orig = data_orig, nimp_cv = nimp_cv,
                        folds = folds, p.crit = p.crit, BW=BW, miceImp = miceImp, ...)
        return(pobjcv)
      }
      
      if(val_method=="cv_MI_RR"){
        if(is_empty(nimp_mice))
          nimp_mice <- pobj$nimp
        pobjcv <- cv_MI_RR(pobj, data_orig = data_orig, nimp_mice = nimp_mice,
                           p.crit = p.crit, BW=BW, folds = folds, miceImp = miceImp, ...)
        return(pobjcv)
      }
    }
    
    if(val_method=="MI_cv_naive"){
      pobjcv <-
        MI_cv_naive(pobj, folds = folds, BW = BW, p.crit = p.crit, cv_naive_appt = cv_naive_appt)
      return(pobjcv)
    }
    
    # Specific for boot_MI and MI_boot
    if(val_method=="boot_MI" | val_method=="MI_boot") {
      if(val_method=="boot_MI"){
        # Part boot_MI
        if(is_empty(data_orig))
          stop("data_orig not defined")
        if(is_empty(nimp_mice))
          nimp_mice <- pobj$nimp
        if (nimp_mice < 2)
          stop("\n", "Number of imputed datasets for MI_boot must be > 1", "\n")
        pobjboot <-
          boot_MI(pobj, data_orig = data_orig, nboot = nboot,
                  nimp_mice = nimp_mice, p.crit = p.crit, direction = direction, 
                  miceImp = miceImp, ...)
        if(p.crit==1)
          message("\n", "p.crit = 1, validation is done without variable selection", "\n")
        return(pobjboot)
      }
      if(val_method=="MI_boot"){
        pobjboot <-
          MI_boot(pobj, p.crit = p.crit, nboot = nboot, direction = direction)
        if(p.crit==1)
          message("\n", "p.crit = 1, validation is done without variable selection", "\n")
        return(pobjboot)
      }
    }
    cal.plot==FALSE
    message("\n", "Calibration plot not made, only when int.val = FALSE", "\n")
  }
}