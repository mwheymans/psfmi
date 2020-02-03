#' Multiparameter pooling methods called by psfmi_mm
#'
#' \code{psfmi_mm_multiparm} Function to pool according to D1, D2 and D3 methods
#'
#' @param data Data frame with stacked multiple imputed datasets.
#'   The original dataset that contains missing values must be excluded from the
#'   dataset. The imputed datasets must be distinguished by an imputation variable,
#'   specified under impvar, and starting by 1 and the clusters should be 
#'   distinguished by a cluster variable, specified under clusvar.
#' @param nimp A numerical scalar. Number of imputed datasets. Default is 5.
#' @param impvar A character vector. Name of the variable that distinguishes the
#'   imputed datasets.
#' @param Outcome Character vector containing the name of the outcome variable.
#' @param P Character vector with the names of the predictor variables.
#'   At least one predictor variable has to be defined.
#' @param p.crit A numerical scalar. P-value selection criterium. A value of 1 
#'   provides the pooled model without selection.
#' @param family Character vector to specify the type of model, "linear" is used to 
#'   call the \code{lmer} function and "binomial" is used to call the \code{glmer}
#'   function of the \code{lme4} package. See details for more information.
#' @param random.eff Character vector to specify the random effects as used by the 
#'   \code{lmer} and \code{glmer} functions of the \code{lme4} package.  
#' @param method A character vector to indicate the pooling method for p-values to pool the
#'   total model or used during predictor selection. This can be "D1", "D2", "D3" or "MPR".
#'   See details for more information.
#' @param print.method logical vector. If TRUE full matrix with p-values of 
#'   all variables according to chosen method (under method) is shown. If FALSE (default) 
#'   p-value for categorical variables according to method are shown and for continuous 
#'   and dichotomous predictors Rubin’s Rules are used.
#'
#' @examples 
#' 
#' \dontrun{
#'  psfmi_mm_multiparm(data=ipdna_md, nimp=5, impvar=".imp", family="linear",
#'  P=c("gender", "bnp", "dbp", "lvef", "bmi_cat"),
#'  random.eff="( 1 | centre)", Outcome="sbp",
#'  p.crit=0.05, method="D1", print.method = FALSE)
#' }
#'
#' @export
psfmi_mm_multiparm <-
  function(data, nimp, impvar, Outcome, P, p.crit, family, random.eff, method, print.method)
{

    pool.p.val <- matrix(0, length(P), 2)
    
    for (j in 1:length(P)) {
      cov.nam0 <- P[-j]
      if (length(P) == 1) {
        cov.nam0 <- "1"
      }
      Y <- c(paste(Outcome, paste("~")))
      form1 <- as.formula(paste(Y, paste(c(P, random.eff), collapse = "+")))
      form0 <- as.formula(paste(Y, paste(c(cov.nam0, random.eff), collapse = "+")))
      
      #coef.fit1 <- se.fit1 <- coef.fit0 <- se.fit0 <- fit.null <- list()
      fit1 <- fit0 <- imp.dt <- list()
      for (i in 1:nimp) {
        imp.dt[[i]] <- data[data[impvar] == i, ]
        fit1[[i]] <- lmer(form1, data = imp.dt[[i]], REML = FALSE)
        fit0[[i]] <- lmer(form0, data = imp.dt[[i]], REML = FALSE)  
        if(family=="binomial"){
        fit1[[i]] <- glmer(form1, data = imp.dt[[i]], family = binomial)
        fit0[[i]] <- glmer(form0, data = imp.dt[[i]], family = binomial)
        } 
      }
      
      out.res1 <- summary(pool(fit1))
      if(family=="binomial"){
      OR <- exp(out.res1[, 1])
      lower.EXP <- exp(out.res1[, 1] - (1.96*out.res1[, 2]))
      upper.EXP <- exp(out.res1[, 1] + (1.96*out.res1[, 2]))
      model.res1 <- cbind(out.res1, OR, lower.EXP, upper.EXP)
      model.res1 <- round(model.res1, 4)
      } else model.res1 <- round(out.res1, 4)
      out.res0 <- summary(mice::pool(fit0))

      tmr <- mitml::testModels(fit1, fit0, method = method)
      
      pvalue <- tmr$test[4]  
      fstat <- tmr$test[1]
      pool.p.val[j, ] <- c(pvalue, fstat)
      
      pool.multiparm <- pool.p.val
      pool.multiparm <- round(data.frame(pool.multiparm), 5)
      row.names(pool.multiparm) <- P
      names(pool.multiparm) <- c("p-values", "F-statistic")
      
    }
    return(pool.multiparm)
}