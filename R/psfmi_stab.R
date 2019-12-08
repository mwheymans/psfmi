#' Function to evaluate predictor and model selection frequency by using bootstrapping
#' over multiple imputed datasets.
#'
#' \code{psfmi_stab} Stability analysis of predictors and prediction models selected with
#'  the \code{psfmi_lr} or \code{psfmi_coxr} functions of the \code{psfmi} package.
#'
#' @param obj An object of class \code{pmimods}, produced by a previous call to \code{psfmi_lr}
#'  or \code{psfmi_coxr}.
#' @param nboot A numerical scalar. Number of bootstrap samples to evaluate the stability. Default is 20.
#'
#' @details The function evaluates predictor selection frequency in stratified bootstrap samples.
#'  The stratification factor is the variable that separates the imputed datasets. It uses as input
#'  an object of class \code{pmimods} as a result of a previous call to the \code{psfmi_lr} or
#'  \code{psfmi_coxr} function.
#'
#'@return A \code{psfmi_stab} object from which the following objects can be extracted: bootstrap
#'  inclusion (selection) frequency of each predictor \code{bif}, total number each predictor is
#'  included in the bootstrap samples as \code{bif_total}, percentage a predictor is selected
#'  in each bootstrap sample as \code{bif_perc} and number of times a prediction model is selected in
#'  the bootstrap samples as \code{model_stab}.
#'
#' @references Heymans MW, van Buuren S. et al. Variable selection under multiple imputation using the bootstrap
#'   in a prognostic study. BMC Med Res Methodol. 2007;13:7-33.
#' @references Eekhout I, van de Wiel MA, Heymans MW. Methods for significance testing of categorical
#'   covariates in logistic regression models after multiple imputation: power and applicability
#'   analysis. BMC Med Res Methodol. 2017;17(1):129.
#' @references Sauerbrei W, Schumacher M. A bootstrap resampling procedure for model building:
#'   application to the Cox regression model. Stat Med. 1992;11:2093–109.
#' @references Royston P, Sauerbrei W (2008) Multivariable model-building – a pragmatic approach to
#'   regression analysis based on fractional polynomials for modelling continuous variables. (2008).
#'   Chapter 8, Model Stability. Wiley, Chichester
#' @references Heinze G, Wallisch C, Dunkler D. Variable selection - A review and
#'  recommendations for the practicing statistician. Biom J. 2018 May;60(3):431-449. 
#'
#' @references http://missingdatasolutions.rbind.io/
#'
#' @examples
#'  pool_lr <- psfmi_lr(data=lbpmilr, nimp=5, impvar="Impnr", Outcome="Chronic",
#'                    predictors=c("Gender", "Smoking",  "JobControl", "JobDemands",
#'                    "Age", "Radiation", "SocialSupport", "Function"),
#'                    cat.predictors = c("Carrying"), p.crit =0.157, method="D1")
#'  pool_lr$RR_Model
#'  pool_lr$multiparm_p
#'
#'  stab_res <- psfmi_stab(pool_lr, nboot=100) # may take a while
#'  stab_res$bif
#'  stab_res$bif_perc
#'  stab_res$model_stab
#'
#' @export
psfmi_stab <- function(obj, nboot=20)
{
  if (class(obj)!="pmimods")
    stop("Object should be of type pmimods")
  call <- match.call()
  nboot <- nboot
  data <- obj$data

  P_select_boot <- list()
  for(i in 1:nboot){
    message("\n", "boot", i, "\n")
    data_boot <- data %>%
      group_by(obj$impvar) %>%
      sample_n(nrow(data), replace = TRUE)
    data_boot <- data.frame(data_boot)[, -ncol(data_boot)]

    if(obj$model_type=="binomial") {
      psfmi_boot <- psfmi_lr(data=data_boot, nimp=obj$nimp, impvar = obj$impvar,
                             Outcome = obj$Outcome, predictors = obj$predictors,
                             p.crit = obj$p.crit, cat.predictors = obj$cat.predictors,
                             spline.predictors = obj$spline.predictors,
                             int.predictors = obj$int.predictors, keep.predictors = obj$keep.predictors,
                             knots = obj$knots, method = obj$method, print.method = obj$print.method)

    }
    if(obj$model_type=="survival") {
      psfmi_boot <- psfmi_coxr(data=data_boot, nimp=obj$nimp, impvar = obj$impvar,
                               time = obj$time, status = obj$status, predictors = obj$predictors,
                               p.crit = obj$p.crit, cat.predictors = obj$cat.predictors,
                               spline.predictors = obj$spline.predictors,
                               int.predictors = obj$int.predictors, keep.predictors = obj$keep.predictors,
                               knots = obj$knots, method = obj$method, print.method = obj$print.method)

    }
    P_select_boot[[i]] <- psfmi_boot$predictors_in[nrow(psfmi_boot$predictors_in), ]
  }

  bif <- data.frame(do.call("rbind", P_select_boot))
  colnames(bif) <- names(P_select_boot[[1]])
  # Group selected models
  bif_pat <- bif %>%
    group_by_all %>%
    count
  # desc order of selected models
  bif_pat_sort <- data.frame(bif_pat %>%
                             arrange(desc(n)) %>%
                             select(1:ncol(bif_pat)))
  bif_pat_perc <- (bif_pat_sort$n / nboot) * 100
  bif_pat_sort <- data.frame(bif_pat_sort, bif_pat_perc)
  colnames(bif_pat_sort) <- c(names(P_select_boot[[1]]), "fr", "perc")

  rownames(bif) <- paste("boot", 1:nboot)

  bif_total <- colSums(bif)
  bif_perc <- (bif_total / nboot) * 100

  stabobj <- list(bif = bif, bif_total = bif_total, bif_perc = bif_perc,
                  model_stab = bif_pat_sort, call = call)
  return(stabobj)
}
