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
#'  pool_performance(data=lbpmilr, nimp=5, impvar="Impnr", 
#'  Outcome = "Chronic", predictors = c("Gender", "Pain", "rcs(Tampascale, 3)", 
#'  "Smoking", "Function", "Radiation", "Age", "factor(Carrying)"), 
#'  cal.plot=TRUE, plot.indiv=FALSE)
#'
#' @export
pool_performance <- function(data, nimp, impvar, Outcome, predictors, 
      cal.plot, plot.indiv, groups_cal=10){

  coef.f <- se.f <- pred.i <- lp.i <- rsq.i <- list()
  roc.f.i <- se.roc.f.i <- se.roc.f.i.logit <- brier.i <- list()
  y.i <- pred.group <- obs.group <- list()

  Y <- c(paste(Outcome, paste("~")))
  fm <- as.formula(paste(Y, paste(predictors, collapse = "+")))
  
  # Pool performance measures over imputed datasets
  for (i in 1:nimp) {

    data_compl <- data[data[impvar] == i, ]
    f <- glm(fm, data = data_compl, family = binomial)

    pred.i[[i]] <- predict(f, type = "response")
    y.i[[i]] <- f$y

    if(cal.plot){
    # Group predicted probabilities
    if(groups_cal ==0) stop("\n", "Number of groups on calibration curve too low", "\n")
    
      group.dec <- cut(pred.i[[i]], quantile(pred.i[[i]],
                                         c(seq(0, 1, 1 / groups_cal))))

    pred.group[[i]] <- tapply(pred.i[[i]], group.dec, mean)
    # Observed probabilities
    obs.group[[i]] <- tapply(y.i[[i]], group.dec, mean)
    }

    lp.i[[i]] <- predict(f)
    coef.f[[i]] <- coef(f)
    se.f[[i]] <- summary(f)[[12]][, 2]

    # ROC/AUC
    roc.f.i[[i]] <- roc(f$y, pred.i[[i]], quiet = TRUE)$auc
    se.roc.f.i[[i]] <- sqrt(pROC::var(roc.f.i[[i]]))
    se.roc.f.i.logit[[i]] <- sqrt(pROC::var(roc.f.i[[i]])) /
      (roc.f.i[[i]]*(1-roc.f.i[[i]]))

    # Nagelkerke R-squared
    f.full <- -2 * logLik(f)
    f.base <- -2 * logLik(update(f, ~ 1))
    n <- f$df.null
    rsq.nagel <- (1 - exp((f.full - f.base)/n))/
      (1 - exp(-f.base/n))
    rsq.i[[i]] <- rsq.nagel

    brier <- function(obs, pred) {
      mean((obs - pred)^2)
    }
    brier.i[[i]] <- brier(f$y, pred.i[[i]])

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
    ((1 + (1/nimp)) * b.roc.logit)
  se.t.roc.logit <- sqrt(tv.roc.logit)

  # Backtransform
  inv.roc <- exp(p.roc.logit) /
    (1 + exp(p.roc.logit))
  inv.roc.u <- exp(p.roc.logit + (1.96*se.t.roc.logit)) /
    (1 + exp(p.roc.logit + (1.96*se.t.roc.logit)))
  inv.roc.l <- exp(p.roc.logit - (1.96*se.t.roc.logit)) /
    (1 + exp(p.roc.logit - (1.96*se.t.roc.logit)))

  roc.m.log <- matrix(c(inv.roc.l, inv.roc, inv.roc.u),
                      1, 3, byrow = T)
  dimnames(roc.m.log) <- list(c("ROC (logit)"),
                              c("95% Low", "ROC", "95% Up"))

  # Pooling R square
  # Fisher z Transformation
  z.rsq <- atanh(unlist(rsq.i))
  z.rsq.p <- mean(z.rsq)

  # within variance
  n <- nrow(data[data[impvar] == 1, ])
  se.z.rsq <- 1/(n-3)
  # between variance
  b.rsq <- var(z.rsq)
  # total variance
  tv.rsq <- se.z.rsq + ((1 + (1/nimp)) * b.rsq)
  se.t.rsq <- sqrt(tv.rsq)
  # inv Fisher z = pooled rsq
  inv.z.rsq.p <- tanh(z.rsq.p)

  # Colmeans of predictors in multiply imputed datasets
  coef_pooled <- colMeans(do.call("rbind", coef.f))

  brier_pool <- mean(unlist(brier.i))

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
  pobjperform <- list(ROC_pooled=roc.m.log, coef_pooled=coef_pooled,
                      R2_pooled=inv.z.rsq.p, Brier_pooled = brier_pool,
                      nimp=nimp)
  # Pooled info in each bootstrap sample
  pobjperform
}
