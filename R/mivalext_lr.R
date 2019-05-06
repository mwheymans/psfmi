#' External Validation of logistic prediction models in MI datasets
#'
#' \code{mivalext_lr} External validation of logistic prediction models
#'
#' @param data.val Data frame or data matrix with stacked multiple imputed validation datasets.
#'  The original dataset that contains missing values must be excluded from the
#'  dataset. The imputed datasets must be distinguished by an imputation variable,
#'  specified under impvar, and starting by 1.
#' @param data.orig A single data frame or data matrix containing the original dataset
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
#'
#' @details The following information of the externally validated model is provided:
#'  pooled ROC curve (median and backtransformed after pooling log transformed
#'  ROC curves), pooled Nagelkerke R-Square value (median and backtransformed
#'  after pooling Fisher transformed values), pooled Hosmer and Lemeshow Test (using
#'  miceadds package), pooled coefficients when model is freely estimated in imputed
#'  datasets and the pooled linear predictor (LP), after the externally validated LP
#'  is estimated in each imputed dataset (provides information about miscalibration
#'  in intercept and slope). When the external validation is very poor,
#'  the R2 fixed can become negative due to the poor fit of the model in
#'  the external dataset (in that case you may report a R2 of zero).
#'
#'@return A \code{mivalext_lr} object from which the following objects 
#'  can be extracted: ROC results as \code{ROC}, R squared results (fixed and calibrated) 
#'  as \code{R2 (fixed)} and \code{R2 (calibr)}, Hosmer and Lemeshow test as \code{HL_test}, 
#'  coefficients pooled as \code{coef_pooled}, linear predictor pooled as \code{LP_pooled ext}, 
#'  and \code{Outcome}, \code{nimp}, \code{impvar}, \code{val.check}, 
#'  \code{g} and \code{coef.check}.
#'
#' @references F. Harrell. Regression Modeling Strategies. With Applications to
#'  Linear Models, Logistic and Ordinal Regression, and Survival Analysis. Springer,
#'  New York, NY, 2015.
#' @references Van Buuren S. (2018). Flexible Imputation of Missing Data. 2nd Edition. Chapman &
#'  Hall/CRC Interdisciplinary Statistics. Boca Raton.
#' @references http://missingdatasolutions.rbind.io/
#' 
#'@examples
#' mivalext_lr(data.val=lbpmilr, nimp=10, impvar="Impnr", Outcome="Chronic",
#' predictors=c("Gender", "factor(Carrying)", "Function", "Tampascale",  "Age"),
#' lp.orig=c(-9.2, -0.34, 0.92, 1.5, 0.5, 0.26, -0.02),
#' cal.plot=TRUE, plot.indiv=TRUE, val.check = TRUE)
#'
#' mivalext_lr(data.val=lbpmilr, nimp=5, impvar="Impnr", Outcome="Chronic",
#' predictors=c("Gender", "factor(Carrying)", "Function", "Tampascale", "Age"),
#' lp.orig=c(-9.2, -0.34, 0.92, 1.1, -0.05, 0.26, -0.02),
#' cal.plot=TRUE, plot.indiv=TRUE, val.check = FALSE)
#'
#' @export
mivalext_lr <-
  function(data.val=NULL, data.orig=NULL, nimp=5, impvar=NULL, Outcome,
   predictors=NULL, lp.orig=NULL, cal.plot=FALSE, plot.indiv=FALSE,
   val.check=FALSE, g=10)
  {
    
    if(is.null(predictors))
      stop("No predictors defined, cannot fit model")
    P <- predictors
    
    # Check data input
    if (!(is.matrix(data.val) | is.data.frame(data.val)))
      stop("Validation dataset should be a matrix or data frame")
    data.val <- data.frame(data.matrix(data.val))
    
    if(!is.null(data.orig)) {
      if (ncol(data.orig) < 2)
        stop("Original data should contain at least two columns")
      if (!(is.matrix(data.orig) | is.data.frame(data.orig)))
        stop("Original dataset should be a matrix or data frame")
      data.orig <- data.frame(data.matrix(data.orig))
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
    if (order(unique(data.val[impvar]))[1] == 0)
      stop("Original dataset should not be included")
    if(is.null(nimp))
      stop("Number of imputed datasets is not defined, use nimp!")
    if (nimp < 2) {
      stop("\n", "Number of imputed datasets must be > 1", "\n\n")
    }
    
    # Determine original (fixed) coefficients
    coef.orig <- lp.orig
    # Determine original coefficients from
    # original dataset
    if(!is.null(data.orig))
    {
      lp.orig <- NULL
      Y <- c(paste(Outcome, paste("~")))
      fm.orig <- as.formula(paste(Y, paste(P, collapse = "+")))
      fit.orig <- glm(fm.orig, x=TRUE, y=TRUE, data=data.orig, family = binomial)
      coef.orig <- coef(fit.orig)
    }
    
    Y <- c(paste(Outcome, paste("~")))
    fm.val <- as.formula(paste(Y, paste(P, collapse = "+")))
    fit.check <- glm(fm.val, x=TRUE, y=TRUE,
                     data=data.val[data.val[impvar] == 1, ], family = binomial)
    coef.check <- names(coef(fit.check))
    # Determine regression formula for correct
    # order of coefficients
    if(val.check==TRUE) {
      #cat("\n", "Order of Predictors, to define lp.orig accordingly", "\n", "\n")
      res.perform <- list("coef.check"=coef.check)
      return(res.perform)
    }
    if (!is.null(lp.orig)){
      if (length(lp.orig) != length(coef.check))
        stop("Number of Predictors not equal to number of coefficients under lp.orig")
    }
    
    rsq.mi.i <- rsq.mi.i.cal <- pred.group <- obs.group <- hl.mi.i <- list()
    roc.f.mi.i <- se.roc.mi.i <- se.roc.mi.i.logit <- list()
    coef.mi.i <- lp.ext.mi <- hl.mi.i <- list()
    
    # Determine performance in each
    # imputed external dataset
    for(i in 1:nimp) {
      data <- data.val[data.val[impvar] == i, ]
      f.ext <- glm(fm.val, data=data, family = binomial)
      X <- model.matrix(f.ext)
      
      lp.ext <- X %*% coef.orig
      f.ext.lp <- glm(f.ext$y ~ lp.ext, family = binomial)
      p.ext <-  c(1/(1+exp(-lp.ext)))
      
      coef.mi.i[[i]] <- coef(f.ext)
      lp.ext.mi[[i]] <- coef(f.ext.lp)
      
      f.ext.stats <- lrm.fit(lp.ext, f.ext$y,
                             initial = c(0, 1), maxit = 1L)
      
      # Nagelkerke R squared
      rsq.mi.i[[i]] <- f.ext.stats$stats["R2"]
      
      # Calibrated R squared
      n <- f.ext.lp$df.null + 1
      k <- f.ext.lp$rank
      logLik1 <- as.numeric(logLik(f.ext.lp))
      f.ext.lp0 <- update(f.ext.lp, . ~ 1)
      logLik0 <- as.numeric(logLik(f.ext.lp0))
      rsq.mi.i.cal[[i]] <- (1 - exp(-2 *
                                      (logLik1 - logLik0)/n)) / (1 - exp(logLik0 * 2/n))
      
      if (cal.plot == TRUE){
        # Group predicted probabilities for calibration curve
        if(length(unique(p.ext))<10) {
          stop("Cannot generate calibration curve, number of groups too small,
        set cal.plot=F")
        } else {
          group.dec <- cut(p.ext, quantile(p.ext,
                                           c(seq(0, 1, 0.1))))
          pred.group[[i]] <- tapply(p.ext, group.dec, mean)
          # Observed probabilities
          obs.group[[i]] <- tapply(f.ext$y, group.dec, mean)
        }
      }
      
      # ROC/AUC
      roc.f.mi.i[[i]] <- roc(f.ext$y, p.ext)$auc
      se.roc.mi.i[[i]] <- sqrt(pROC::var(roc.f.mi.i[[i]]))
      se.roc.mi.i.logit[[i]] <- sqrt(pROC::var(roc.f.mi.i[[i]])) /
        (roc.f.mi.i[[i]]*(1-roc.f.mi.i[[i]]))
      
      # Hosmer and Lemeshow Chi square value
      if(g<4){
        stop("For Hosmer and Lemeshow test, number of groups must be > 3")
      } else {
        hl.mi.i[[i]] <- hoslem.test(f.ext$y, p.ext, g=g)[[1]]
      }
    }
    
    coef.pool <- round(colMeans(do.call("rbind", coef.mi.i)), 5)
    lp.pool <- round(colMeans(do.call("rbind", lp.ext.mi)), 5)
    
    # ROC/AUC
    # RR on logit transformation ROC curve and SE
    est.roc.logit <- log(unlist(roc.f.mi.i)/
                           (1-unlist(roc.f.mi.i)))
    se.roc.logit <- unlist(se.roc.mi.i.logit)
    
    # Pooling
    p.roc.logit <- mean(est.roc.logit)
    # within variance
    p.se.roc.logit <- mean(se.roc.logit)
    # between variance
    b.roc.logit <- var(est.roc.logit)
    # total variance
    tv.roc.logit <- p.se.roc.logit +
      ((1 + (1/nimp)) * b.roc.logit)
    se.t.roc.logit <- sqrt(tv.roc.logit)
    
    # Backtransform
    inv.roc <- exp(p.roc.logit) /
      (1 + exp(p.roc.logit))
    inv.roc.u <- exp(p.roc.logit + (1.96*se.t.roc.logit)) /
      (1 + exp(p.roc.logit + (1.96*se.t.roc.logit)))
    inv.roc.l <- exp(p.roc.logit - (1.96*se.t.roc.logit)) /
      (1 + exp(p.roc.logit - (1.96*se.t.roc.logit)))
    
    roc.m.log <- round(matrix(c(inv.roc.l, inv.roc, inv.roc.u),
                              1, 3, byrow = TRUE), 5)
    dimnames(roc.m.log) <- list(c("ROC (logit)"),
                                c("95% Low", "ROC", "95% Up"))
    
    # Median and IQR ROC
    roc.med.iqr <- round(summary(unlist(roc.f.mi.i))[-c(1, 4, 6)], 5)
    
    roc.res <- list("ROC (logit)"=roc.m.log,
                    "ROC (median)"=roc.med.iqr)
    
    #### Pooling R square (uncalibrated)
    # Fisher z Transformation
    z.rsq <- atanh(unlist(rsq.mi.i))
    z.rsq.p <- mean(z.rsq)
    
    # within variance
    n <- nrow(data.val[data.val[, impvar] == 1, ])
    se.z.rsq <- 1/(n-3)
    # between variance
    b.rsq <- var(z.rsq)
    # total variance
    tv.rsq <- se.z.rsq + ((1 + (1/nimp)) * b.rsq)
    se.t.rsq <- sqrt(tv.rsq)
    # inv Fisher z = pooled rsq
    inv.z.rsq.p <- round(tanh(z.rsq.p), 5)
    
    # Median and IQR R square
    rsq.med.iqr <- round(summary(unlist(rsq.mi.i))[-c(1,4,6)], 5)
    
    res.rsq <- list("Fisher Z (fixed)"=inv.z.rsq.p,
                    "Median (fixed)"=rsq.med.iqr)
    
    #### Pooling R square (calibrated)
    # Fisher z Transformation
    z.rsq.cal <- atanh(unlist(rsq.mi.i.cal))
    z.rsq.p.cal <- mean(z.rsq.cal)
    
    # within variance
    n <- nrow(data.val[data.val[, impvar] == 1, ])
    se.z.rsq.cal <- 1/(n-3)
    # between variance
    b.rsq.cal <- var(z.rsq.cal)
    # total variance
    tv.rsq.cal <- se.z.rsq.cal + ((1 + (1/nimp)) * b.rsq.cal)
    se.t.rsq.cal <- sqrt(tv.rsq.cal)
    # inv Fisher z = pooled rsq
    inv.z.rsq.p.cal <- round(tanh(z.rsq.p.cal), 5)
    
    # Median and IQR R square
    rsq.med.iqr.cal <- round(summary(unlist(rsq.mi.i.cal))[-c(1,4,6)], 5)
    
    res.rsq.cal <- list("Fisher Z (calibrated)"=inv.z.rsq.p.cal,
                        "Median (calibrated)"=rsq.med.iqr.cal)
    
    # H&L test
    res.hl <- round(miceadds::micombine.chisquare(unlist(hl.mi.i),
                                                  g-2, display = F), 5)
    
    message("\n", "Pooled performance measures over m = ",
        nimp, " imputed external validation datasets correctly estimated", "\n\n")
    
    res.perform <- list("ROC"=roc.res, "R2 (fixed)"=res.rsq,
      "R2 (calibr)"=res.rsq.cal, "HL test"=res.hl, "coef_pooled"=coef.pool,
      "LP_pooled ext"=lp.pool, "nimp"=nimp, "impvar"=impvar,
      "Outcome"=Outcome, "val.check"=val.check, "g"=g, "coef.check"=coef.check)
    
    if(cal.plot==TRUE) {
      ID.mi <- rep(1:nimp, each=10)
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