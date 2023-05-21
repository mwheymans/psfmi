#' Function to evaluate bootstrap predictor and model stability.
#'
#' \code{stab_single} Stability analysis of predictors and prediction models selected with
#'  the \code{glm_bw}.
#'
#' @param pobj An object of class \code{smods} (single models), produced by a previous call to
#'  \code{glm_bw}.
#' @param nboot A numerical scalar. Number of bootstrap samples to evaluate the stability. Default is 20.
#' @param p.crit A numerical scalar. Used as P-value selection criterium during bootstrap model selection.
#' @param start_model If TRUE the bootstrap evaluation takes place from the start model of object pobj, if
#'  FALSE the final model is used for the evaluation.
#'
#' @details The function evaluates predictor selection frequency in bootstrap samples.
#'  It uses as input an object of class \code{smods} as a result of a
#'  previous call to the \code{glm_bw}.
#'
#'@return A \code{psfmi_stab} object from which the following objects can be extracted: bootstrap
#'  inclusion (selection) frequency of each predictor \code{bif}, total number each predictor is
#'  included in the bootstrap samples as \code{bif_total}, percentage a predictor is selected
#'  in each bootstrap sample as \code{bif_perc} and number of times a prediction model is selected in
#'  the bootstrap samples as \code{model_stab}.
#'
#' @references Heymans MW, van Buuren S. et al. Variable selection under multiple imputation using the bootstrap
#'   in a prognostic study. BMC Med Res Methodol. 2007;13:7-33.
#' @references Sauerbrei W, Schumacher M. A bootstrap resampling procedure for model building:
#'   application to the Cox regression model. Stat Med. 1992;11:2093–109.
#' @references Royston P, Sauerbrei W (2008) Multivariable model-building – a pragmatic approach to
#'   regression analysis based on fractional polynomials for modelling continuous variables. (2008).
#'   Chapter 8, Model Stability. Wiley, Chichester.
#' @references Heinze G, Wallisch C, Dunkler D. Variable selection - A review and
#'  recommendations for the practicing statistician. Biom J. 2018;60(3):431-449.
#'
#' @references http://missingdatasolutions.rbind.io/
#'
#' @examples
#'  model_lr <- glm_bw(formula = Radiation ~ Pain + factor(Satisfaction) + 
#'    rcs(Tampascale,3) + Age + Duration + JobControl + JobDemands + SocialSupport, 
#'    data=lbpmilr_dev, p.crit = 0.05)
#'
#' \dontrun{
#'  stab_res <- stab_single(model_lr, start_model = TRUE, nboot=20, p.crit=0.05)
#'  stab_res$bif
#'  stab_res$bif_perc
#'  stab_res$model_stab
#'}
#'
#' @export
stab_single <- function(pobj, nboot=20, p.crit = 0.05, start_model = TRUE)
{
  if(!inherits(pobj, "smods"))
    stop("\n", "Object should be of type smods", "\n")
  if(p.crit==1)
    stop("\n", "To determine Model Stability p.crit must be < 1", "\n\n")
  if(sum(is.na(pobj$data))!=0)
    stop("\n", "Bootstrap samples can not contain missing values", "\n\n")
  if(start_model == FALSE & is_empty(pobj$predictors_final))
    stop( "\n", "Final model is empty. Set start_model=TRUE to select from the initial model", "\n\n")
  
  call <- match.call()
  nboot <- nboot
  data <- pobj$data
  boot_seq <- as.list(1:nboot)
  
  boot_data <- bootstraps(data, times = nboot)
  boot_pred_pat <- mapply(function(x, y) {
    message("\n", "Boot ", y)
    x <- as.data.frame(x)
    
    if(start_model){
      psfmi_boot <- glm_bw(formula = pobj$formula_initial, data = x,
                              p.crit = p.crit, keep.predictors = pobj$keep.predictors,
                              model_type = pobj$model_type)
    } else {
      psfmi_boot <- glm_bw(formula = pobj$formula_final, data = x,
                              p.crit = p.crit, keep.predictors = pobj$keep.predictors,
                              model_type = pobj$model_type)
    }
    
    boot_predictors_in <- ifelse(psfmi_boot$predictors_out[nrow(psfmi_boot$predictors_out), ], 0, 1)
    return(boot_predictors_in)
  }, x = boot_data$splits, y=boot_seq, SIMPLIFY = FALSE)
  
  bif <- data.frame(do.call("rbind", boot_pred_pat))
  
  if(!start_model)
    colnames(bif) <- pobj$predictors_final
  if(start_model)
    colnames(bif) <- pobj$predictors_initial
  
  # Group selected models
  bif_pat <- bif %>%
    group_by_all %>%
    count
  
  # desc order of selected models
  bif_pat_sort <- 
    data.frame(bif_pat[order(bif_pat$n, decreasing = TRUE), ])
  bif_pat_perc <- 
    round((bif_pat_sort$n / nboot) * 100, 0)
  bif_pat_sort <- 
    data.frame(bif_pat_sort, bif_pat_perc)
  
  if(!start_model) {
    colnames(bif_pat_sort) <- c(pobj$predictors_final, "freq", "bif_pat_perc")
  } else {
    colnames(bif_pat_sort) <- c(pobj$predictors_initial, "freq", "bif_pat_perc")
  }
  rownames(bif) <- paste("boot", 1:nboot)
  
  bif_total <- colSums(bif)
  bif_perc <- round((bif_total / nboot) * 100, 3)
  
  stabobj <- list(bif = bif, bif_total = bif_total, bif_perc = bif_perc,
                  model_stab = bif_pat_sort)
  return(stabobj)
}