#' Pools the Likelihood Ratio tests across Multiply Imputed datasets ( method D4) 
#'
#' \code{pool_D4} The D4 statistic to combine the likelihood ratio tests (LRT)
#'  across Multiply Imputed datasets according method D4.
#'  
#' @param data Data frame with stacked multiple imputed datasets.
#'   The original dataset that contains missing values must be excluded from the
#'   dataset. The imputed datasets must be distinguished by an imputation variable,
#'   specified under impvar, and starting by 1.
#' @param nimp A numerical scalar. Number of imputed datasets. Default is 5.
#' @param impvar A character vector. Name of the variable that distinguishes the
#' imputed datasets. 
#' @param fm0 the null model.
#' @param fm1 the (nested) model to compare. Must be larger than the null model.
#' @param robust if TRUE a robust LRT is used (algorithm 1 in Chan and Meng), otherwise 
#'  algorithm 2 is used. 
#' @param model_type if TRUE (default) a logistic regression model is fitted, otherwise
#'  a linear regression model is used
#'                                                                                                                                                                                            
#' @return The D4 statistic, the numerator, df1 and denominator, df2 degrees of freedom 
#'  for the F-test.
#' 
#' @references Chan, K. W., & Meng, X.-L. (2019). Multiple improvements of multiple 
#'  imputation likelihood ratio tests. ArXiv:1711.08822 [Math, Stat]. https://arxiv.org/abs/1711.08822
#' @references Grund, Simon, Oliver Lüdtke, and Alexander Robitzsch. 2021. “Pooling Methods for 
#'  Likelihood Ratio Tests in Multiply Imputed Data Sets.” PsyArXiv. January 29. doi:10.31234/osf.io/d459g.
#'   
#' @author Martijn Heymans, 2021
#'
#' @examples
#'
#' fm0 <- Chronic ~ BMI + factor(Carrying) + 
#'   Satisfaction + SocialSupport + Smoking
#' fm1 <- Chronic ~ BMI + factor(Carrying) + 
#'   Satisfaction +  SocialSupport + Smoking +
#'   Radiation
#'
#' psfmi::pool_D4(data=lbpmilr, nimp=10, impvar="Impnr",
#'                fm0=fm0, fm1=fm1, robust = TRUE)
#'    
#' @export 
pool_D4 <- function(data,
                    nimp,
                    impvar,
                    fm0,
                    fm1,
                    robust=TRUE,
                    model_type="binomial")
{
  m <- nimp
  data <-
    filter(data, data[[impvar]] <= nimp)
  
  if(model_type=="binomial") family="binomial"
  if(model_type=="linear") family="gaussian"
  
  ll0 <- ll1 <- list()
  for (i in 1:nimp) {
    imp.dt <- data[data[impvar] == i, ]
    fit0 <- glm(fm0, data = imp.dt, family = family)
    fit1 <- glm(fm1, data = imp.dt, family = family)
    ll1[[i]] <- logLik(fit1)
    ll0[[i]] <- logLik(fit0)
  }
  
  # Average over imputed datasets
  d_L_bar <-
    mean(-2*(mapply(function(x, y)
      x-y, x=ll0, y=ll1)))
  
  d_L_bar
  
  # Log likelihood in stacked data
  ll0_S <-
    logLik(glm(fm0,
               family = family, data=data)) / nimp
  
  ll1_S <-
    logLik(glm(fm1,
               family = family, data=data)) / nimp
  
  # pooled LR in stacked data
  d_S_hat <- -2 * (ll0_S - ll1_S)
  
  df0 <-
    attr(ll0[[1]], "df")
  df1 <-
    attr(ll1[[1]], "df")
  h <- df1
  k <- df1-df0
  
  if(robust){ # algorithm 1
    deltabar <- 2 * mean(unlist(ll1)) # input 3.7
    deltahat <- 2 * ll1_S # input 3.7
    r <- (m + 1) / (h * (m - 1)) * (deltabar - deltahat)
  } else { # Algorithm 2
    rs <- (m + 1) / (k * (m - 1)) * (d_L_bar - d_S_hat)
    r <- max(0, rs)
  }
  
  # D4
  D <- d_S_hat / (k * (1 + r))
  
  if(robust){
    df2 <- ((1 + r)/r)^2 * h * (m-1)
  } else {
    df2 <- ((1 + r)/r)^2 * k * (m-1)
    
  }
  
  pval <- pf(D, k, df2, lower.tail = FALSE)
  objd4 <- list(F = D[1], df1 = k[1],
                df2 = df2[1], pval = pval[1], riv = r[1])
  return(objd4)
}