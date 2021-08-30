#' Pooling performance measures across multiply imputed datasets
#'
#' \code{pool_performance} Pooling performance measures for logistic
#'  regression models.
#'
#' @param data Data frame with stacked multiple imputed datasets.
#'   The original dataset that contains missing values must be excluded 
#'   from the dataset.
#' @param formula A formula object to specify the model as normally used by glm.
#' @param nimp A numerical scalar. Number of imputed datasets. Default is 5.
#' @param impvar A character vector. Name of the variable that distinguishes 
#'   the imputed datasets.
#' @param model_type If "binomial" (default), performance measures are calculated
#'  for logistic regression models, if "survival" for Cox regression models.      
#' @param cal.plot If TRUE a calibration plot is generated. Default is TRUE. 
#'  model_type must be "binomial".
#' @param plot.method If "mean" one calibration plot is generated, first taking the 
#'   mean of the linear predictor across the multiply imputed datasets (default), if 
#'   "individual" the calibration plot of each imputed dataset is plotted, 
#'   if "overlay" calibration plots from each imputed datasets are plotted in one figure. 
#' @param plot.indiv This argument is deprecated; please use plot.method instead.  
#' @param groups_cal A numerical scalar. Number of groups used on the calibration plot and. 
#'  for the Hosmer and Lemeshow test. Default is 10. If the range of predicted probabilities. 
#'  is low, less than 10 groups can be chosen, but not < 3. 
#'
#' @examples
#'  perf <- pool_performance(data=lbpmilr, nimp=5, impvar="Impnr", 
#'  formula = Chronic ~ Gender + Pain + Tampascale + 
#'  Smoking + Function + Radiation + Age + factor(Carrying), 
#'  cal.plot=TRUE, plot.method="mean", 
#'  groups_cal=10, model_type="binomial")
#'  
#'  perf$ROC_pooled
#'  perf$R2_pooled
#'  
#' @export
pool_performance <- function(data, 
                             formula, 
                             nimp, 
                             impvar, 
                             model_type="binomial",
                             cal.plot=TRUE, 
                             plot.method="mean", 
                             plot.indiv,
                             groups_cal=10)
 {
  
  coef_f <- se_f <- pred.group <- obs.group <- list()
  
  if (!missing(plot.indiv)) {
    warning("argument plot.indiv is deprecated; please use plot.method instead.",
            call. = FALSE)
    plot.method <- plot.indiv
  }
  if(model_type=="survival" & cal.plot)
    stop("Calibration plots only possible for logistic regression models, change model_type")
  
  if(is_empty(formula))
    stop("\n", "Model not specified in formula object")
  
  fm <-
    formula
  
  if(!any(grepl("Surv", fm)) & model_type=="survival")
    stop("When model_type is survival, outcome in formula must be a Surv object")
  if(any(grepl("Surv", fm)) & model_type=="binomial")
    stop("When outcome is a Surv object, model_type must be of type survival")
  
  if(model_type=="binomial"){
    perf_stats <-
      matrix(NA, nimp, 5)
  } else { perf_stats <-
    matrix(NA, nimp, 3) }
  
  lp_mi <- matrix(NA, nrow(data[data[impvar] == 1, ]), nimp)
  
  # Pool performance measures over imputed datasets
  for (i in 1:nimp) {
    
    data_compl <- data[data[impvar] == i, ]
    
    if(model_type=="binomial"){
      f <-
        glm(fm, data = data_compl, family = binomial)
      pred <-
        predict(f, type = "response")
      
      lp_mi[, i] <- predict(f)
      
      if(cal.plot){
        # Group predicted probabilities
        if(groups_cal <=3) stop("\n", "Number of groups on calibration curve must be > 3", "\n")
        group.dec <- cut(pred, quantile(pred,
                                        c(seq(0, 1, 1 / groups_cal))))
        # Predicted probabilities
        pred.group[[i]] <-
          tapply(pred, group.dec, mean)
        # Observed probabilities
        obs.group[[i]] <-
          tapply(f$y, group.dec, mean)
      }
      
      coef_f[[i]] <-
        coef(f)
      se_f[[i]] <-
        summary(f)[[12]][, 2]
      
      # ROC/AUC
      roc.i <-
        pROC::roc(f$y, pred, quiet = TRUE)$auc
      se.roc.i <-
        sqrt(pROC::var(roc.i))
     
      # Nagelkerke R-squared
      rsq.i <-
        rsq_nagel(f)
      sc_brier.i <-
        scaled_brier(f$y, pred)

      hoslem_chi <- hoslem_test(f$y, pred, g=groups_cal)[[1]]
      
      perf_stats[i, ] <-
        c(roc.i, se.roc.i, rsq.i, sc_brier.i, hoslem_chi)
    }
    
    if(model_type=="survival"){
      f <-
        coxph(fm, data = data_compl)
      
      coef_f[[i]] <-
        coef(f)
      se_f[[i]] <-
        summary(f)[[7]][, 3]
      
      # concordance
      c_stats <-
        f$concordance[6]
      se_c_stats <-
        f$concordance[7]
      
      # Nagelkerke R-squared
      rsq.i <-
        rsq_surv(f)
      perf_stats[i, ] <-
        c(c_stats, se_c_stats, rsq.i)
    }
  }
  # End pooling performance measures in multiply imputed data
  
  if(model_type=="binomial"){
    hltest_pooled <- 
      pool_D2(perf_stats[, 5], v=groups_cal-2)
    sc_brier_pool <-
      mean(perf_stats[, 4])
  }
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
        tapply(f$y, group.dec, mean)
      mean.cal.plot <- data.frame("Obs"=obs.group,
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
  # Pooling Performance measures and coefficients
  if(model_type=="binomial"){
    pobjperform <- list(ROC_pooled=roc_res, coef_pooled=coef_pooled,
                        R2_pooled=rsq.n, Brier_Scaled_pooled = sc_brier_pool,
                        nimp=nimp, HLtest_pooled=hltest_pooled,
                        model_type = model_type)
  } else {
    pobjperform <- list(concordance_pooled=roc_res, coef_pooled=coef_pooled,
                        R2_pooled=rsq.n, nimp=nimp,
                        model_type = model_type)
  }
  # Pooled info in each bootstrap sample
  return(pobjperform)
}