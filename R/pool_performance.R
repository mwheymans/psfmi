#' Pooling performance measures over multiply imputed datasets
#'
#' \code{pool_performance} Pooling performance measures
#'
#' @param data Data frame with stacked multiple imputed datasets.
#'   The original dataset that contains missing values must be excluded 
#'   from the dataset.
#' @param nimp A numerical scalar. Number of imputed datasets. Default is 5.
#' @param impvar A character vector. Name of the variable that distinguishes 
#'   the imputed datasets.
#' @param Outcome Character vector containing the name of the outcome variable.
#' @param predictors Character vector with the names of the predictor 
#'   variables as used in the formula part of an glm object.
#' @param cal.plot If TRUE a calibration plot is generated. Default is FALSE. 
#'   Can be used in combination with int_val = FALSE.
#' @param plot.indiv If TRUE calibration plots for each separate imputed dataset 
#'   are generated, otherwise all calibration plots are plotted in one figure.       
#' @param groups_cal A numerical scalar. Number of groups used on the calibration plot. 
#'  Default is 10. If the range of predicted probabilities is low, less than 10 groups 
#'  can be chosen. 
#'
#'@examples 
#'  perf <- pool_performance(data=lbpmilr, nimp=5, impvar="Impnr", 
#'  Outcome = "Chronic", predictors = c("Gender", "Pain", "rcs(Tampascale, 3)", 
#'  "Smoking", "Function", "Radiation", "Age", "factor(Carrying)"), 
#'  cal.plot=TRUE, plot.indiv=FALSE)
#'  
#'  perf$ROC_pooled
#'  
#' @export
pool_performance <- function(data, nimp, impvar, Outcome, predictors,
                             cal.plot, plot.indiv, groups_cal=10){
  
  coef_f <- se_f <- pred.group <- obs.group <- list()
  
  Y <-
    c(paste(Outcome, paste("~")))
  fm <-
    as.formula(paste(Y, paste(predictors, collapse = "+")))
  
  perf_stats <-
    matrix(NA, nimp, 4)
  
  # Pool performance measures over imputed datasets
  for (i in 1:nimp) {
    
    data_compl <- data[data[impvar] == i, ]
    f <-
      glm(fm, data = data_compl, family = binomial)
    P <-
      predict(f, type = "response")
    
    if(cal.plot){
      # Group predicted probabilities
      if(groups_cal ==0) stop("\n", "Number of groups on calibration curve too low", "\n")
      group.dec <- cut(P, quantile(P,
                                   c(seq(0, 1, 1 / groups_cal))))
      # Predicted probabilities
      pred.group[[i]] <-
        tapply(P, group.dec, mean)
      # Observed probabilities
      obs.group[[i]] <-
        tapply(f$y, group.dec, mean)
    }
    
    #LP[[i]] <- predict(f)
    coef_f[[i]] <-
      coef(f)
    se_f[[i]] <-
      summary(f)[[12]][, 2]
    
    # ROC/AUC
    roc.i <-
      pROC::roc(f$y, P, quiet = TRUE)$auc
    se.roc.i <-
      sqrt(pROC::var(roc.i))
    
    # Nagelkerke R-squared
    rsq.i <-
      rsq_nagel(f)
    sc_brier.i <-
      scaled_brier(f$y, P)
    
    perf_stats[i, ] <-
      c(roc.i, se.roc.i, rsq.i, sc_brier.i)
  }
  # End pooling performance measures in multiply imputed data
  
  # ROC/AUC
  roc_res <-
    pool_auc(perf_stats[, 1], perf_stats[, 2],
             nimp = nimp, log_auc = TRUE)
  
  # Pooling R square
  # Fisher z Transformation
  rsq.n <-
    tanh(mean(atanh(perf_stats[, 3])))
  
  # Colmeans of predictors in multiply imputed datasets
  coef_pooled <-
    colMeans(do.call("rbind", coef_f))
  
  sc_brier_pool <-
    mean(perf_stats[, 4])
  
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
    if(plot.indiv){
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
  # Pooling Performance measures and coefficients
  pobjperform <- list(ROC_pooled=roc_res, coef_pooled=coef_pooled,
                      R2_pooled=rsq.n, Brier_Scaled_pooled = sc_brier_pool,
                      nimp=nimp)
  # Pooled info in each bootstrap sample
  pobjperform
}