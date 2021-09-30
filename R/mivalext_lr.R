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
#' @param formula A formula object to specify the model as normally used by glm.
#' @param lp.orig Numeric vector of the original coefficient values that are
#'  externally validated.
#' @param cal.plot If TRUE a calibration plot is generated. Default is FALSE.
#' @param plot.indiv This argument is deprecated; please use plot.method instead.  
#' @param val.check logical vector. If TRUE the names of the predictors of the LP
#'  are provided and can be used as information for the order of the coefficient
#'  values as input for lp.orig. If FALSE (default) validation procedure is executed
#'  with coefficient values fitted in the order as used under lp.orig.
#' @param g A numerical scalar. Number of groups for the Hosmer and
#'  Lemeshow test. Default is 10.
#' @param groups_cal A numerical scalar. Number of groups used on the calibration plot. 
#'  Default is 10. If the range of predicted probabilities is low, less than 10 groups 
#'  can be chosen.
#' @param plot.method If "mean" one calibration plot is generated, first taking the 
#'   mean of the linear predictor values across the multiply imputed datasets (default), if 
#'   "individual" the calibration plot in each imputed dataset is plotted, 
#'   if "overlay" calibration plots from each imputed datasets are plotted in one figure. 
#'
#' @details The following information of the externally validated model is provided:
#'  \code{ROC} pooled ROC curve (back transformed after pooling log transformed
#'  ROC curves), \code{R2} pooled Nagelkerke R-Square value (back transformed after 
#'  pooling Fisher transformed values), \code{HLtest} pooled Hosmer and Lemeshow 
#'  Test (using function \code{pool_D2}), \code{coef_pooled} pooled coefficients 
#'  when model is freely estimated in imputed datasets and \code{LP_pooled_ext} 
#'  the pooled linear predictor (LP), after the externally validated LP is estimated in 
#'  each imputed dataset (provides information about miscalibration in intercept and slope). 
#'  In addition information is provided about \code{nimp}, \code{impvar}, \code{formula},
#'  \code{val_ckeck}, \code{g} and \code{coef_check}. When the external validation is 
#'  very poor, the R2 can become negative due to the poor fit of the model in
#'  the external dataset (in that case you may report a R2 of zero).
#'
#'@return A \code{mivalext_lr} object from which the following objects 
#'  can be extracted: ROC results as \code{ROC}, R squared results as \code{R2}, 
#'  Hosmer and Lemeshow test as \code{HL_test}, coefficients pooled as 
#'  \code{coef_pooled}, linear predictor pooled as \code{LP_pooled ext}, 
#'  and \code{formula}, \code{nimp}, \code{impvar}, \code{val.check}, 
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
#' mivalext_lr(data.val=lbpmilr, nimp=5, impvar="Impnr", 
#'   formula = Chronic ~ Gender + factor(Carrying)  + Function + 
#'   Tampascale + Age, lp.orig=c(-10, -0.35, 1.00, 1.00, -0.04, 0.26, -0.01),
#'   cal.plot=TRUE, val.check = FALSE)
#'
#' @export
mivalext_lr <-
  function(data.val=NULL,
           data.orig=NULL,
           nimp=5,
           impvar=NULL,
           formula=NULL,
           lp.orig=NULL,
           cal.plot=FALSE,
           plot.indiv,
           val.check=FALSE,
           g=10,
           groups_cal=10,
           plot.method="mean")
  {
    
    if (!missing(plot.indiv)) {
      warning("argument plot.indiv is deprecated; please use plot.method instead.",
              call. = FALSE)
      plot.method <- plot.indiv
    }
    if(is_empty(formula))
      stop("\n", "Model not specified in formula object")
    form <- terms(formula)
    Outcome <-
      as.character(attr(form, "variables")[[2]])
    Y <-
      c(paste(Outcome, paste("~")))
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
    if(!is.null(data.orig))  {
      lp.orig <-
        NULL
      
      fm.orig <-
        formula
      
      fit.orig <-
        glm(fm.orig, x=TRUE, y=TRUE,
            data=data.orig, family = binomial)
      coef.orig <-
        coef(fit.orig)
    }
    
    fm.val <-
      formula
    fit.check <-
      glm(fm.val, x=TRUE, y=TRUE,
          data=data.val[data.val[impvar] == 1, ],
          family = binomial)
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
    
    LLlogistic <-
      function(formula, data, coefs) {
        logistic <- function(mu) exp(mu)/(1 + exp(mu))
        Xb <- model.matrix(formula, data) %*% coefs
        y <- model.frame(formula, data)[1][, 1]
        p <- logistic(Xb)
        y <- (y - min(y))/(max(y) - min(y))
        term1 <- term2 <- rep(0, length(y))
        term1[y != 0] <- y[y != 0] * log(y[y != 0]/p[y != 0])
        term2[y == 0] <- (1 - y[y == 0]) * log((1 - y[y == 0])/(1 - p[y == 0]))
        return(-(2 * sum(term1 + term2)))
      }
    
    pred.group <- obs.group <- coef_extern <- list()
    
    stats_ext <-
      matrix(NA, nrow = nimp, ncol = 12)
    
    n <-
      nrow(data.val[data.val[impvar] == 1, ])
    
    lp_mi <-
      matrix(NA, n, nimp)
    
    # Determine performance in each
    # imputed external dataset
    for(i in 1:nimp) {
      data <-
        data.val[data.val[impvar] == i, ]
      #n <-
      #  nrow(data)
      f.ext <-
        glm(fm.val, data=data, family = binomial)
      X <-
        model.matrix(f.ext)
      
      lp.ext <-
        X %*% coef.orig
      lp_mi[, i] <-
        lp.ext
      f.ext.lp <-
        glm(f.ext$y ~ lp.ext, family = binomial)
      p.ext <-
        c(1/(1+exp(-lp.ext)))
      
      coef_extern[[i]] <-
        coef(f.ext)
      slope_extern <-
        coef(f.ext.lp)
      slope_extern_se <-
        sqrt(diag(vcov(f.ext.lp)))
      
      lp_offset_int <-
        glm(f.ext$y ~ offset(lp.ext), family = binomial)
      lp_offset_int_se <-
        sqrt(diag(vcov(lp_offset_int)))
      lp_offset_slope <-
        glm(f.ext$y ~ lp.ext + offset(lp.ext), family = binomial)
      lp_offset_slope_se <-
        sqrt(diag(vcov(lp_offset_slope)))
      
      fit_full <-
        -1*LLlogistic(fm.val, data = data,
                      coef.orig)
      
      fit_null <-
        -2*logLik(glm(as.formula(paste(Y, paste(1))),
                      data = data, family=binomial))
      
      
      # Nagelkerke R squared
      rsq.nagel <-
        c((1 - exp((fit_full - fit_null)/n))/
            (1 - exp(-fit_null/n)))[1]
      #if(rsq.nagel<0) rsq.nagel <- 0.00000001
      
      if(cal.plot){
        # Group predicted probabilities
        if(groups_cal <4) stop("\n", "Number of groups on calibration curve must be > 3", "\n")
        group.dec <- cut(p.ext, quantile(p.ext,
                                         c(seq(0, 1, 1 / groups_cal))))
        # Predicted probabilities
        pred.group[[i]] <-
          tapply(p.ext, group.dec, mean)
        # Observed probabilities
        obs.group[[i]] <-
          tapply(f.ext$y, group.dec, mean)
      }
      
      # ROC/AUC
      auc.m <-
        roc(f.ext$y, p.ext, quiet = TRUE)$auc
      se.roc.m <-
        sqrt(pROC::var(auc.m))
      
      # Hosmer and Lemeshow Chi square value
      if(g<4){
        stop("For Hosmer and Lemeshow test, number of groups must be > 3")
      } else {
        hl.m <- hoslem_test(f.ext$y, p.ext, g=g)$chisq
      }
      
      stats_ext[i,] <-  c(slope_extern, slope_extern_se, rsq.nagel, auc.m, se.roc.m, hl.m,
                          coef(lp_offset_int), lp_offset_int_se, coef(lp_offset_slope)[2],
                          lp_offset_slope_se[2])
    }
    
    stats_ext <-
      data.frame(intercept=stats_ext[, 1], intercept_se=stats_ext[, 3], slope=stats_ext[, 2],
                 slope_se=stats_ext[, 4], rsq.mi=stats_ext[, 5],
                 auc.mi=stats_ext[, 6], auc.mi.se=stats_ext[, 7],
                 hl.mi=stats_ext[, 8], offset_intercept=stats_ext[, 9],
                 offset_intercept_se=stats_ext[, 10], offset_slope=stats_ext[, 11],
                 offset_slope_se=stats_ext[, 12])
    
    coef.pool <-
      round(colMeans(do.call("rbind", coef_extern)), 5)
    
    source("pool_RR.R")
    
    pooled_int <-
      pool_RR(stats_ext$intercept, stats_ext$intercept_se, n=n, k=1)
    
    pooled_slope <-
      pool_RR(stats_ext$slope, stats_ext$slope_se, n=n, k=1)
    pooled_offset_int <-
      pool_RR(stats_ext$offset_intercept, stats_ext$offset_intercept_se, n=n, k=1)
    pooled_offset_slope <-
      pool_RR(stats_ext$offset_slope, stats_ext$offset_slope_se, n=n, k=1)
    res_pool_lp <-
      rbind(pooled_int, pooled_slope, pooled_offset_int, pooled_offset_slope)
    
    # ROC/AUC
    # RR on logit transformation ROC curve and SE
    auc_RR <-
      pool_auc(stats_ext$auc.mi, stats_ext$auc.mi.se,
               nimp = nimp, log_auc = TRUE)
    
    # Fisher z Transformation
    z.rsq <-
      atanh(stats_ext$rsq.mi)
    z.rsq.p <-
      mean(z.rsq)
    # inv Fisher z = pooled rsq
    inv.z.rsq.p <- round(tanh(z.rsq.p), 5)
    
    # H&L test
    res.hl <-
      round(pool_D2(dw=stats_ext[, 6], v=g-2), 5)
    
    message("\n", "Pooled performance measures over m = ",
            nimp, " imputed external validation datasets
            correctly estimated", "\n\n")
    
    res.perform <-
      list(Calibrate_LP=res_pool_lp, coef_pooled=coef.pool,
           ROC=auc_RR, R2=inv.z.rsq.p, HLtest=res.hl,
           nimp=nimp, impvar=impvar, formula=formula,
           val_check=val.check, g=g,
           coef_check=coef.check,
           groups_cal=groups_cal)
    if(cal.plot){
      ID.mi <- rep(1:nimp, each=groups_cal)
      myX <- scale_x_continuous(limits = c(-0.1, 1.1),
                                breaks=seq(0,1,0.1),
                                name = "Predicted Probabilities")
      myY <- scale_y_continuous(limits = c(-0.1, 1.1),
                                breaks=seq(0,1,0.1),
                                name = "Observed Probabilities")
      data.cal.plot <- data.frame(ID.mi, "Obs"=unlist(obs.group),
                                  "Pred"=unlist(pred.group))
      theme_set(theme_bw())
      if(plot.method=="mean")
      {
        pred_mean <-
          1 / ( 1 + exp(-rowMeans(lp_mi)))
        # Group predicted probabilities
        group.dec <- cut(pred_mean, quantile(pred_mean,
                                             c(seq(0, 1, 1 / groups_cal))))
        # Predicted probabilities
        pred.group <-
          tapply(pred_mean, group.dec, mean)
        
        # Observed probabilities
        obs.group <-
          do.call("rbind", obs.group)
        
        obs.group <- ifelse(obs.group==0, 0.000001, obs.group)
        obs.group <- ifelse(obs.group==1, 0.999999, obs.group)
        
        log_mean_obs_group <-
          rowMeans(apply(obs.group,1, FUN = function (x) (log(x / (1-x)))))
        obs.group_exp <- exp(log_mean_obs_group) / (1 + exp(log_mean_obs_group))
        
        #obs.group <-
        #  colMeans(do.call("rbind", obs.group))
        mean.cal.plot <- data.frame("Obs"=obs.group_exp,
                                    "Pred"=pred.group)
        g1 <- ggplot(data = mean.cal.plot, aes_string(x = "Pred",
                                                      y = "Obs")) + geom_point() +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()) + myX + myY
        g2 <- g1 + stat_smooth(method = "lm", se = FALSE,
                               formula = y ~ splines::bs(x, 3))
        g3 <- g2 + geom_abline(slope=1, intercept=0, linetype="dashed")
        print(g3)
      }
      if(plot.method=="individual"){
        # Calibration plot in each imputed dataset
        g1 <- ggplot(data = data.cal.plot, aes_string(x = "Pred", y = "Obs",
                                                      group = "ID.mi")) + geom_point() +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())
        g2 <- g1 + stat_smooth(method = "lm", se = FALSE,
                               formula = y ~ splines::bs(x, 3)) +
          facet_wrap(~ ID.mi) + myX + myY
        g3 <- g2 + geom_abline(slope=1, intercept=0, linetype="dashed")
        print(g3)
      }
      if(plot.method=="overlay") {
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