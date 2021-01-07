#' External Validation of logistic prediction models in 
#'  multiply imputed datasets
#'
#' \code{mivalext_lr} External validation of logistic prediction models
#'
#' @param data.val Data frame with stacked multiply imputed validation datasets.
#'  The original dataset that contains missing values must be excluded from the
#'  dataset. The imputed datasets must be distinguished by an imputation variable,
#'  specified under impvar, and starting by 1.
#' @param data.orig A single data frame containing the original dataset
#'  that was used to develop the model. Used to estimate the original regression
#'  coefficients in case lp.orig is not provided.
#' @param nimp A numerical scalar. Number of imputed datasets. Default is 5.
#' @param impvar A character vector. Name of the variable that distinguishes the
#'  imputed datasets.
#' @param Outcome Character vector containing the name of the outcome variable.
#' @param predictors Character vector with the names of the predictor variables
#'  of the model that is validated.
#' @param lp.orig Numeric vector of the original coefficient values that are
#'  externally validated.
#' @param cal.plot If TRUE a calibration plot is generated. Default is FALSE.
#' @param plot.indiv If TRUE calibration plots of each imputed dataset are
#'  generated. Default is FALSE.
#' @param val.check logical vector. If TRUE the names of the predictors of the LP
#'  are provided and can be used as information for the order of the coefficient
#'  values as input for lp.orig. If FALSE (default) validation procedure is executed
#'  with coefficient values fitted in the order as used under lp.orig.
#' @param g A numerical scalar. Number of groups for the Hosmer and
#'  Lemeshow test. Default is 10.
#' @param groups_cal A numerical scalar. Number of groups used on the calibration plot. 
#'  Default is 10. If the range of predicted probabilities is low, less than 10 groups 
#'  can be chosen.
#'
#' @details The following information of the externally validated model is provided:
#'  \code{ROC} pooled ROC curve (median and back transformed after pooling log transformed
#'  ROC curves), \code{R2_fixed} and \code{R2_calibr} pooled Nagelkerke R-Square value 
#'  (median and back transformed after pooling Fisher transformed values), \code{HLtest} 
#'  pooled Hosmer and Lemeshow Test (using miceadds package), \code{coef_pooled} pooled 
#'  coefficients when model is freely estimated in imputed datasets and \code{LP_pooled_ext} 
#'  the pooled linear predictor (LP), after the externally validated LP is estimated in 
#'  each imputed dataset (provides information about miscalibration in intercept and slope). 
#'  In addition information is provided about \code{nimp}, \code{impvar}, \code{Outcome},
#'  \code{val_ckeck}, \code{g} and \code{coef_check}. When the external validation is 
#'  very poor, the R2 fixed can become negative due to the poor fit of the model in
#'  the external dataset (in that case you may report a R2 of zero).
#'
#'@return A \code{mivalext_lr} object from which the following objects 
#'  can be extracted: ROC results as \code{ROC}, R squared results (fixed and calibrated) 
#'  as \code{R2 (fixed)} and \code{R2 (calibr)}, Hosmer and Lemeshow test as \code{HL_test}, 
#'  coefficients pooled as \code{coef_pooled}, linear predictor pooled as \code{LP_pooled ext}, 
#'  and \code{Outcome}, \code{nimp}, \code{impvar}, \code{val.check}, 
#'  \code{g}, \code{coef.check} and \code{groups_cal}.
#'
#' @references F. Harrell. Regression Modeling Strategies. With Applications to
#'  Linear Models, Logistic and Ordinal Regression, and Survival Analysis. 2nd Edition.
#'  Springer, New York, NY, 2015.
#' @references Van Buuren S. (2018). Flexible Imputation of Missing Data. 2nd Edition. Chapman &
#'  Hall/CRC Interdisciplinary Statistics. Boca Raton.
#' @references http://missingdatasolutions.rbind.io/
#' 
#'@examples
#'
#' mivalext_lr(data.val=lbpmilr, nimp=5, impvar="Impnr", Outcome="Chronic",
#' predictors=c("Gender", "factor(Carrying)", "Function", "Tampascale", "Age"),
#' lp.orig=c(-10, -0.35, 1.00, 1.00, -0.04, 0.26, -0.01),
#' cal.plot=TRUE, plot.indiv=TRUE, val.check = FALSE)
#'
#' @export
mivalext_lr <-
  function(data.val=NULL, 
           data.orig=NULL, 
           nimp=5, 
           impvar=NULL, 
           Outcome,
           predictors=NULL, 
           lp.orig=NULL, 
           cal.plot=FALSE, 
           plot.indiv=FALSE,
           val.check=FALSE, 
           g=10, 
           groups_cal=10)
 {
    
    if(is.null(predictors))
      stop("No predictors defined, cannot fit model")
    P <-
      predictors
    
    # Check data input
    if (!(is.data.frame(data.val)))
      stop("Data should be a data frame")
    data.val <-
      data.frame(as_tibble(data.val))
    data.val <-
      mutate_if(data.val, is.factor, ~ as.numeric(as.character(.x)))
    
    if(!is.null(data.orig)) {
      if (ncol(data.orig) < 2)
        stop("Original data should contain at least two columns")
      if (!(is.data.frame(data.orig)))
        stop("Original dataset should be a data frame")
      data.orig <- data.frame(as_tibble(data.orig))
      data.orig <- mutate_if(data.orig, is.factor, ~ as.numeric(as.character(.x)))
      if (!is.null(lp.orig))
        stop("Not needed to define lp.orig, coefficient values
      are derived from original dataset as defined under data.orig")
    }
  
    if(is.null(data.orig) & is.null(lp.orig))
      stop("data.orig and lp.orig not defined, no model to validate")
    if (ncol(data.val) < 2)
      stop("Validation data should contain at least two columns")
    if(is.null(impvar))
      stop("Imputation variable is not defined")
    if (sort(unique(data.val[,impvar]))[1] == 0)
      stop("Original dataset should not be included")
    if(is.null(nimp))
      stop("Number of imputed datasets is not defined, use nimp!")
    if (nimp < 2) {
      stop("\n", "Number of imputed datasets must be > 1", "\n\n")
    }
    
    # Determine original (fixed) coefficients
    coef.orig <-
      lp.orig
    # Determine original coefficients from
    # original dataset
    if(!is.null(data.orig))
    {
      lp.orig <-
        NULL
      Y <-
        c(paste(Outcome, paste("~")))
      fm.orig <-
        as.formula(paste(Y, paste(P, collapse = "+")))
      fit.orig <-
        glm(fm.orig, x=TRUE, y=TRUE, data=data.orig, family = binomial)
      coef.orig <-
        coef(fit.orig)
    }
    
    Y <-
      c(paste(Outcome, paste("~")))
    fm.val <-
      as.formula(paste(Y, paste(P, collapse = "+")))
    fit.check <-
      glm(fm.val, x=TRUE, y=TRUE,
          data=data.val[data.val[impvar] == 1, ], family = binomial)
    coef.check <-
      names(coef(fit.check))
    # Determine regression formula for correct
    # order of coefficients
    if(val.check==TRUE) {
      res.perform <-
        list("coef.check"=coef.check)
      return(res.perform)
    }
    if (!is.null(lp.orig)){
      if (length(lp.orig) != length(coef.check))
        stop("Number of Predictors not equal to number of coefficients under lp.orig")
    }
    
    pred.group <- obs.group <- coef.mi <- list()
    
    stats_ext <-
      matrix(NA, nrow = nimp, ncol = 7)
    # Determine performance in each
    # imputed external dataset
    for(i in 1:nimp) {
      data <-
        data.val[data.val[impvar] == i, ]
      f.ext <-
        glm(fm.val, data=data, family = binomial)
      X <-
        model.matrix(f.ext)
      
      lp.ext <-
        X %*% coef.orig
      f.ext.lp <-
        glm(f.ext$y ~ lp.ext, family = binomial)
      p.ext <-
        c(1/(1+exp(-lp.ext)))
      
      coef.mi[[i]] <-
        coef(f.ext)
      lp_ext <-
        coef(f.ext.lp)
      
      f.ext.stats <-
        lrm.fit(lp.ext, f.ext$y, initial = c(0, 1), maxit = 1L)
      
      # Nagelkerke R squared
      rsq.mi <-
        f.ext.stats$stats["R2"]
      
      # Calibrated R squared
      n <-
        f.ext.lp$df.null + 1
      k <-
        f.ext.lp$rank
      logLik1 <-
        as.numeric(logLik(f.ext.lp))
      f.ext.lp0 <-
        update(f.ext.lp, . ~ 1)
      logLik0 <-
        as.numeric(logLik(f.ext.lp0))
      rsq.mi.cal <-
        (1 - exp(-2 * (logLik1 - logLik0)/n)) / (1 - exp(logLik0 * 2/n))
      
      if (cal.plot){
        # Group predicted probabilities for calibration curve
        if(groups_cal ==0) stop("\n", "Number of groups on calibration curve too low", "\n")
        group.dec <- cut(p.ext, quantile(p.ext,
                                         c(seq(0, 1, 1 / groups_cal))))
        pred.group[[i]] <- tapply(p.ext, group.dec, mean)
        # Observed probabilities
        obs.group[[i]] <- tapply(f.ext$y, group.dec, mean)
      }
      
      # ROC/AUC
      auc.mi <-
        roc(f.ext$y, p.ext, quiet = TRUE)$auc
      se.roc.mi <-
        sqrt(pROC::var(auc.mi))

      # Hosmer and Lemeshow Chi square value
      if(g<4){
        stop("For Hosmer and Lemeshow test, number of groups must be > 3")
      } else {
        hl.mi <- hoslem.test(f.ext$y, p.ext, g=g)[[1]]
      }
      
      stats_ext[i,] <-  c(lp_ext, rsq.mi, rsq.mi.cal, auc.mi, se.roc.mi, hl.mi)
    }
    
    stats_ext <- data.frame(intercept=stats_ext[, 1], slope=stats_ext[, 2], rsq.mi=stats_ext[, 3],
                            rsq.mi.cal=stats_ext[, 4], auc.mi=stats_ext[, 5], auc.mi.se=stats_ext[, 6],
                            hl.mi=stats_ext[, 7])
    
    coef.pool <- round(colMeans(do.call("rbind", coef.mi)), 5)
    lp.pool <- round(colMeans(stats_ext[, c("intercept", "slope")]), 5)
    
    # ROC/AUC
    # RR on logit transformation ROC curve and SE
    auc_RR <- pool_auc(stats_ext$auc.mi, stats_ext$auc.mi.se,
                       nimp = nimp, log_auc = TRUE)
    
    # Median and IQR ROC
    roc.med.iqr <- round(summary(stats_ext$auc.mi)[-c(1, 4, 6)], 5)
    
    roc.res <- list("ROC (logit)"=auc_RR,
                    "ROC (median)"=roc.med.iqr)
    
    # Pooling R square (uncalibrated)
    # Fisher z Transformation
    z.rsq <- atanh(stats_ext$rsq.mi)
    z.rsq.p <- mean(z.rsq)
    
    # inv Fisher z = pooled rsq
    inv.z.rsq.p <- round(tanh(z.rsq.p), 5)
    
    # Median and IQR R square
    rsq.med.iqr <- round(summary(stats_ext$rsq.mi)[-c(1,4,6)], 5)
    
    res.rsq <- list("Fisher Z (fixed)"=inv.z.rsq.p,
                    "Median (fixed)"=rsq.med.iqr)
    
    # Pooling R square (calibrated)
    # Fisher z Transformation
    z.rsq.cal <- atanh(stats_ext$rsq.mi.cal)
    z.rsq.p.cal <- mean(z.rsq.cal)

    # inv Fisher z = pooled rsq
    inv.z.rsq.p.cal <- round(tanh(z.rsq.p.cal), 5)
    
    # Median and IQR R square
    rsq.med.iqr.cal <- round(summary(stats_ext$rsq.mi.cal)[-c(1,4,6)], 5)
    
    res.rsq.cal <- list("Fisher Z (calibrated)"=inv.z.rsq.p.cal,
                        "Median (calibrated)"=rsq.med.iqr.cal)
    
    # H&L test
    res.hl <- round(miceadds::micombine.chisquare(stats_ext$hl.mi,
                                                  g-2, display = F), 5)
    
    message("\n", "Pooled performance measures over m = ",
            nimp, " imputed external validation datasets correctly estimated", "\n\n")
    
    res.perform <- list(ROC=roc.res, R2_fixed=res.rsq,
                        R2_calibr=res.rsq.cal, HLtest=res.hl, coef_pooled=coef.pool,
                        LP_pooled_ext=lp.pool, nimp=nimp, impvar=impvar,
                        Outcome=Outcome, val_check=val.check, g=g, coef_check=coef.check,
                        groups_cal=groups_cal)
    
    if(cal.plot==TRUE) {
      ID.mi <- rep(1:nimp, each = groups_cal)
      myX <- scale_x_continuous(limits = c(-0.1, 1.1),
                                breaks=seq(0,1,0.1),
                                name = "Predicted Probabilities")
      myY <- scale_y_continuous(limits = c(-0.1, 1.1),
                                breaks=seq(0,1,0.1),
                                name = "Observed Probabilities")
      data.cal.plot <- data.frame(ID.mi, "Obs"=unlist(obs.group),
                                  "Pred"=unlist(pred.group))
      theme_set(theme_bw())
      if(plot.indiv==TRUE){
        # Calibration plot in each imputed dataset
        g1 <- ggplot(data = data.cal.plot, aes_string(x = "Pred", y = "Obs",
                                                      group = "ID.mi")) + geom_point() + theme(panel.grid.major = element_blank(),
                                                                                               panel.grid.minor = element_blank())
        g2 <- g1 + stat_smooth(method = "lm", se = FALSE,
                               formula = y ~ splines::bs(x, 3)) +
          facet_wrap(~ ID.mi) + myX + myY
        g3 <- g2 + geom_abline(slope=1, intercept=0, linetype="dashed")
        print(g3)
      } else {
        # Overlaying Calibration plots
        g1 <- ggplot(data = data.cal.plot, aes_string(x = "Pred",
                                                      y = "Obs", group = "ID.mi")) + geom_point() +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()) + myX + myY
        g2 <- g1 + stat_smooth(method = "lm", se = FALSE,
                               formula = y ~ splines::bs(x, 3))
        g3 <- g2 + geom_abline(slope=1, intercept=0, linetype="dashed")
        print(g3)
      }
    }
    return(res.perform)
}