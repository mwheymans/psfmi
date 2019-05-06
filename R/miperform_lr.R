#' Evaluate performance of logistic regression models over MI datasets
#'
#' \code{miperform_lr} Evaluate Performance of logistic regression models
#'
#' @param data Data frame or data matrix with stacked multiple imputed datasets.
#'   The original dataset that contains missing values must be excluded from the
#'   dataset. The imputed datasets must be distinguished by an imputation variable,
#'   specified under impvar.
#' @param nimp A numerical scalar. Number of imputed datasets. Default is 5.
#' @param impvar A character vector. Name of the variable that distinguishes the
#'   imputed datasets.
#' @param Outcome Character vector containing the name of the outcome variable.
#' @param predictors Character vector with the names of the predictor variables.
#'   At least one predictor variable must be defined.
#' @param cat.predictors A single string or a vector of strings to define the
#'   categorical variables. Default is NULL categorical predictors.
#' @param int.predictors A single string or a vector of strings with the names of
#'   the variables that form an interaction pair, separated by a “:” symbol.
#' @param cal.plot If TRUE a calibration plot is generated. Default is FALSE.
#' @param plot.indiv If TRUE calibration plots of each imputed dataset are
#' generated.
#' @param int.val If TRUE performance measures are reported as a result
#'  of internal validation in each imputed datasets. This is a wrapper function
#'  of Frank Harrell´s validate function as part of the rms package.
#' @param method "boot" is the default setting to generate bootstrap corrected
#'  performance measures.
#' @param B The number of bootstrap resamples, default is 250.
#' @param bw If TRUE backward selection is applied during bootstrap
#'  internal validation. Default is FALSE. Backward selection is done
#'  using the fastbw function of the rms package.
#' @param rule Set at "p" for backward selection using the p-value as
#'  criterium when bw=TRUE.
#' @param type Set at "individual" for backward selection of
#'  individual predictors when bw=TRUE.
#' @param p.val P-value criterium for backward selection when bw=TRUE.
#' @param force A vector of integers to define the variables that are forced
#'  in the model during backward selection. The integer value matches
#'  the order of the variable in the model (starting with the intercept).
#'  
#'@return A \code{miperform_lr} object from which the following objects 
#'  can be extracted: ROC results as \code{ROC}, R squared results as \code{R2}, 
#'  Hosmer and Lemeshow test as \code{HL_test}, linear predictor pooled as \code{LP_pooled}, 
#'  performance after internal validation as \code{Int_val_pooled}, 
#'  and \code{Outcome}, \code{nimp}, \code{impvar}, \code{predictors}, 
#'  \code{cat.predictors}, \code{int.predictors}, \code{int.val}.
#'  
#' @references Marshall A, Altman DG, Holder RL, Royston P. Combining estimates of
#'  interest in prognostic modelling studies after multiple imputation: current
#'  practice and guidelines. BMC Med Res Methodol. 2009;9:57.
#' @references F. Harrell. Regression Modeling Strategies. With Applications to
#'  Linear Models, Logistic and Ordinal Regression, and Survival Analysis. Springer,
#'  New York, NY, 2015.
#' @references Van Buuren S. (2018). Flexible Imputation of Missing Data. 2nd Edition. Chapman &
#'  Hall/CRC Interdisciplinary Statistics. Boca Raton.
#' @references Harel, O. (2009). The estimation of R2 and adjusted R2 in
#'  incomplete data sets using multiple imputation. Journal of Applied Statistics,
#'  36(10), 1109-1118
#' @references http://missingdatasolutions.rbind.io/
#' 
#'@examples
#' miperform_lr(data=lbpmilr, nimp=5, impvar="Impnr",
#' Outcome=c("Chronic"), predictors=c("Gender", "Pain",
#' "Tampascale","Smoking","Function", "Radiation", "Age"),
#' cat.predictors=c("Carrying", "Satisfaction"),
#' int.predictors=c("Carrying:Smoking", "Gender:Smoking"),
#' cal.plot=TRUE, plot.indiv = FALSE)
#'
#' @export
miperform_lr <-
  function(data, nimp=5, impvar=NULL, Outcome,
   predictors=NULL, cat.predictors=NULL, int.predictors=NULL,
   cal.plot=FALSE, plot.indiv=FALSE, int.val=FALSE,
   method="boot", B=250, bw=FALSE, rule="p", type="individual",
   p.val=0.05, force=NULL)
  {
    P <- predictors
    cat.P <- cat.predictors
    int.P <- int.predictors
    P.check <-c(P, cat.P)
    
    # Check data input
    if (!(is.matrix(data) | is.data.frame(data)))
      stop("Data should be a matrix or data frame")
    data <- data.frame(data.matrix(data))
    if ((nvar <- ncol(data)) < 2)
      stop("Data should contain at least two columns")
    if(is.null(impvar))
      stop("Imputation variable is not defined")
    if (order(unique(data[impvar]))[1] == 0)
      stop("Original dataset should not be included")
    if(is.null(nimp))
      stop("Number of imputed datasets is not defined, use nimp!")
    if (nimp < 2) {
      stop("\n", "Number of imputed datasets must be > 1", "\n\n")
    }
    if (!is.null(cat.P)) {
      if(any(cat.P%in%P)){
        cat.P.double <- cat.P[cat.P%in%P]
        stop("\n", "Categorical variable(s) -", cat.P.double,
             "- also defined as Predictor", "\n\n")
      }
    }
    if(any(duplicated(P))){
      stop("\n", "Predictor(s) - ", c(P[duplicated(P)]),
           " - defined more than once", "\n\n")
    }
    # Check if al variables are available in dataset
    if(any(!P.check %in% names(data))) {
      P.mis <- P.check[!P.check %in% names(data)]
      stop("\n", "Predictor(s) - ", P.mis,
           "- not available in dataset", "\n\n")
    }
    if(!is.null(int.P)) {
      int.P.check <- lapply(int.P[grep(":", int.P)],
                            function(x) { unlist(strsplit(x, split=":")) })
      int.P.check <- unique(unlist(int.P.check))
      if(any(!int.P.check %in% P.check))
        stop("\n", "Not all interaction terms defined as
          Predictor or Categorical Predictor", "\n\n")
    }
    # First predictors, second cetegorical
    # predictors and last interactions!
    P <- c(P, cat.P, int.P)
    if (is.null(P))
      stop("\n", "No predictors defined, cannot fit model", "\n\n")
    
    Y <- c(paste(Outcome, paste("~")))
    # Define formula for validate function
    fm.int <- as.formula(paste(Y, paste(P, collapse = "+")))
    
    if (!is.null(cat.P)) {
      if(length(cat.P)==1){
        P <- gsub(cat.P,
                  replacement=paste0("factor(", cat.P, ")"), P)
        data[[cat.P]] <- factor(data[[cat.P]])
      } else {
        for(i in 1:length(cat.P)) {
          P <- gsub(cat.P[i],
                    replacement=paste0("factor(", cat.P[i], ")"), P)
          data[[cat.P[[i]]]] <- factor(data[[cat.P[[i]]]])
        }
      }
    }
    
    levels.cat.P <- lapply(cat.P, function(x) {
      nr.levels.cat.P <- length(table(data[data[impvar] == 1, ][, x]))
      if (nr.levels.cat.P < 3) {
        stop("\n", "Categorical variable(s) only 2 levels,
          do not define as categorical", "\n\n")
      }
    })
    
    # Start  loop for backward selection over imputed datasets
    y.i <- coef.f <- se.f <- pred.i <- lp.i <- rsq.i <- hl.i <- list()
    roc.f.i <- se.roc.f.i <- se.roc.f.i.logit <- list()
    int.val.i <- auc.iv.i <- rsq.iv.i <- slope.iv.i <- list()
    pred.group <- obs.group <- list()
    
    fm <- as.formula(paste(Y, paste(P, collapse = "+")))
    
    # Calculate pooled LP
    data.imp <- split(data[data[,impvar] %in% c(1:nimp), ],
                      data[data[,impvar] %in% c(1:nimp), impvar])
    p.model <- lapply(data.imp, function(x) {
      coef(glm(fm, family=binomial, data=x))
    })
    lp.pool <- colMeans(do.call("rbind", p.model))
    
    # Start loop for Rubin's Rules
    for (i in 1:nimp) {
      f <- glm(fm, data = data[data[impvar] == i, ],
               family = binomial)
      y.i[[i]] <- f$y
      pred.i[[i]] <- predict(f, type = "response")
      
      # Group predicted probabilities
      group.dec <- cut(pred.i[[i]], quantile(pred.i[[i]],
                                             c(seq(0, 1, 0.1))))
      pred.group[[i]] <- tapply(pred.i[[i]], group.dec, mean)
      # Observed probabilities
      obs.group[[i]] <- tapply(y.i[[i]], group.dec, mean)
      
      lp.i[[i]] <- predict(f)
      coef.f[[i]] <- coef(f)
      se.f[[i]] <- summary(f)[[12]][, 2]
      
      # Nagelkerke R-squared
      f.full <- -2 * logLik(f)
      f.base <- -2 * logLik(update(f, ~ 1))
      n <- f$df.null
      rsq.nagel <- (1 - exp((f.full - f.base)/n))/
        (1 - exp(-f.base/n))
      rsq.i[[i]] <- rsq.nagel
      
      # ROC/AUC
      roc.f.i[[i]] <- roc(f$y, pred.i[[i]])$auc
      se.roc.f.i[[i]] <- sqrt(pROC::var(roc.f.i[[i]]))
      se.roc.f.i.logit[[i]] <- sqrt(pROC::var(roc.f.i[[i]])) /
        (roc.f.i[[i]]*(1-roc.f.i[[i]]))
      
      # Hosmer and Lemeshow test
      hl.i[[i]] <- hoslem.test(f$y, pred.i[[i]])[[1]]
      
      #### Internal bootstrap validation
      if(int.val==T) {
        f <- rms::lrm(fm.int, data = data[data[impvar] == i, ], x=T, y=T)
        int.val.i[[i]] <- rms::validate(f, method=method, B=B, bw=bw,
         rule=rule, type=type, sls=p.val, force=force, estimates=F)
        auc.iv.i[[i]] <- (int.val.i[[i]][45]+1)/2
        rsq.iv.i[[i]] <- int.val.i[[i]][46]
        slope.iv.i[[i]] <- int.val.i[[i]][48]
      }
    }
    
    # ROC/AUC
    # RR on logit transformation ROC curve and SE
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
    
    # Median and IQR ROC
    roc.med.iqr <- summary(unlist(roc.f.i))[-c(1, 4, 6)]
    
    roc.res <- list("ROC (logit)"=roc.m.log,
                    "ROC (median)"=roc.med.iqr)
    
    #### Pooling R square
    # Fisher z Transformation
    z.rsq <- atanh(unlist(rsq.i))
    z.rsq.p <- mean(z.rsq)
    
    # within variance
    n <- nrow(data[data[, impvar] == 1, ])
    se.z.rsq <- 1/(n-3)
    # between variance
    b.rsq <- var(z.rsq)
    # total variance
    tv.rsq <- se.z.rsq + ((1 + (1/nimp)) * b.rsq)
    se.t.rsq <- sqrt(tv.rsq)
    # inv Fisher z = pooled rsq
    inv.z.rsq.p <- tanh(z.rsq.p)
    
    # Median and IQR R square
    rsq.med.iqr <- summary(unlist(rsq.i))[-c(1,4,6)]
    
    res.rsq <- list("R2 (Fisher Z)"=inv.z.rsq.p,
                    "Median R2"=rsq.med.iqr)
    
    # H&L test
    res.hl <- miceadds::micombine.chisquare(unlist(hl.i),
      8, display = F)
    
    message("\n", "Pooled performance measures over m = ",
            nimp, " imputed datasets correctly estimated", "\n")
    res.perform <- list("ROC"=roc.res, "R2"=res.rsq,
      "HL_test"=res.hl, "LP_pooled"=lp.pool,
      "nimp"=nimp, "impvar"=impvar, "Outcome"=Outcome,
      "predictors"=predictors, "cat.predictors"=cat.predictors,
      "int.predictors"=int.predictors, "int.val"=int.val)
    
    if(int.val==T){
      # Result Internal validation
      int.val.roc <- summary(unlist(auc.iv.i))[-c(1,6)]
      int.val.rsq <- summary(unlist(rsq.iv.i))[-c(1,6)]
      int.val.slope <- summary(unlist(slope.iv.i))[-c(1,6)]
      int.val.res <- round(rbind(int.val.roc, int.val.rsq,
                                 int.val.slope), 4)
      row.names(int.val.res) <- c("ROC", "R-Square", "Slope")
      int.val.res <- list("Bootstrap corrected"=int.val.res)
      message("\n", "Internal validation using ",B,
              " Bootstrap samples correctly estimated", "\n\n")
      
      res.perform <- list("ROC"=roc.res, "R2"=res.rsq,
        "HL_test"=res.hl, "LP_pooled"=lp.pool, "Int_val_pooled"=int.val.res, 
        "nimp"=nimp, "impvar"=impvar, "Outcome"=Outcome,
        "predictors"=predictors, "cat.predictors"=cat.predictors,
        "int.predictors"=int.predictors, "int.val"=int.val)
    }
    
    if(cal.plot==T) {
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
      if(plot.indiv==T){
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
