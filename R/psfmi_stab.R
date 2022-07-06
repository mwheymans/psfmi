#' Function to evaluate bootstrap predictor and model stability in multiply imputed datasets.
#'
#' \code{psfmi_stab} Stability analysis of predictors and prediction models selected with
#'  the \code{psfmi_lr}, \code{psfmi_coxr} or \code{psfmi_mm} functions of the \code{psfmi} package.
#'
#' @param pobj An object of class \code{pmods} (pooled models), produced by a previous call to
#'  \code{psfmi_lr}, \code{psfmi_coxr} or \code{psfmi_mm}.
#' @param boot_method A single string to define the bootstrap method. Use "single" after a call to
#'  \code{psfmi_lr} and \code{psfmi_coxr} and "cluster" after a call to \code{psfmi_mm}.
#' @param nboot A numerical scalar. Number of bootstrap samples to evaluate the stability. Default is 20.
#' @param p.crit A numerical scalar. Used as P-value selection criterium during bootstrap model selection.
#' @param start_model If TRUE the bootstrap evaluation takes place from the start model of object pobj, if
#'  FALSE the final model is used for the evaluation.
#' @param direction The direction of predictor selection, "BW" for backward selection and "FW"
#'   for forward selection.
#'#'
#' @details The function evaluates predictor selection frequency in stratified or cluster bootstrap samples.
#'  The stratification factor is the variable that separates the imputed datasets. The same bootstrap cases
#'  are drawn in each bootstrap sample. It uses as input an object of class \code{pmods} as a result of a
#'  previous call to the \code{psfmi_lr}, \code{psfmi_coxr} or \code{psfmi_mm} functions.
#'  In combination with the \code{psfmi_mm} function a cluster bootstrap method is used where bootstrapping
#'  is used on the level of the clusters only (and not also within the clusters).
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
#'  recommendations for the practicing statistician. Biom J. 2018;60(3):431-449.
#'
#' @references http://missingdatasolutions.rbind.io/
#'
#' @section Vignettes:
#'  https://mwheymans.github.io/psfmi/articles/psfmi_StabilityAnalysis.html
#'
#' @examples
#'  pool_lr <- psfmi_coxr(formula = Surv(Time, Status) ~ Pain + factor(Satisfaction) + 
#'    rcs(Tampascale,3) + Radiation + Radiation*factor(Satisfaction) + Age + Duration + 
#'    Previous + Radiation*rcs(Tampascale, 3), data=lbpmicox, p.crit = 0.157, direction="FW",
#'    nimp=5, impvar="Impnr", keep.predictors = NULL, method="D1")
#'
#'  pool_lr$RR_Model
#'  pool_lr$multiparm
#'
#' \dontrun{
#'  stab_res <- psfmi_stab(pool_lr, direction="FW", start_model = TRUE,
#'      boot_method = "single", nboot=20, p.crit=0.05)
#'  stab_res$bif
#'  stab_res$bif_perc
#'  stab_res$model_stab
#'}
#'
#' @export
psfmi_stab <- function(pobj, boot_method=NULL, nboot=20,
                       p.crit = 0.05, start_model = TRUE, direction = NULL)
{
  
  if(!inherits(pobj, "pmods"))
    stop("\n", "Object should be of type pmods", "\n")
  if(is.null(boot_method))
    stop("\n", "boot_method is not defined, choose single or cluster", "\n")
  if(p.crit==1)
    stop("\n", "To determine Model Stability p.crit must be < 1", "\n\n")
  
  if(boot_method=="single" & pobj$model_type=="binomial" | pobj$model_type=="survival"){
    if(is_empty(direction))
      stop( "\n", "Specify FW or BW for forward or backward predictor selection", "\n")
    if(start_model == FALSE & is_empty(pobj$predictors_final))
      stop( "\n", "Final model is empty. You cannot determine the stability of an empty model", "\n\n")
    if(pobj$model_type=="binomial"){
      if(pobj$direction=="FW" & start_model == FALSE){
        Y <-
          c(paste(pobj$Outcome, paste("~")))
        pobj$formula_final <-
          as.formula(paste(Y, paste(pobj$predictors_final, collapse = "+")))
      }
    }
    if(pobj$model_type=="survival"){
      if(pobj$direction=="FW" & start_model == FALSE){
        Y <-
          c(paste0("Surv(", pobj$time, ",", pobj$status, ")~"))
        pobj$formula_final <-
          as.formula(paste(Y, paste(pobj$predictors_final, collapse = "+")))
      }
    }
  }
  call <- match.call()
  nboot <- nboot
  data <- pobj$data
  boot_seq <- as.list(1:nboot)
  
  if (inherits(pobj, "pmods") & boot_method == "single"){
    
    boot_data <- bootstraps(data, strata = pobj$impvar, times = nboot)
    boot_pred_pat <- mapply(function(x, y) {
      message("\n", "Boot ", y)
      x <- as.data.frame(x)
      
      if(pobj$model_type=="binomial") {
        
        if(start_model){
          psfmi_boot <- psfmi_lr(formula = pobj$formula_initial, data = x, nimp=pobj$nimp, impvar = pobj$impvar,
                                 p.crit = p.crit, keep.predictors = pobj$keep.predictors,
                                 method = pobj$method, direction = direction)
          # print(psfmi_boot)
          
        } else {
          psfmi_boot <- psfmi_lr(formula = pobj$formula_final, data = x, nimp=pobj$nimp, impvar = pobj$impvar,
                                 p.crit = p.crit, keep.predictors = pobj$keep.predictors,
                                 method = pobj$method, direction = direction)
        }
      }
      if(pobj$model_type=="survival") {
        
        if(start_model){
          psfmi_boot <- psfmi_coxr(formula = pobj$formula_initial, data = x, nimp=pobj$nimp, impvar = pobj$impvar,
                                   p.crit = p.crit, keep.predictors = pobj$keep.predictors,
                                   method = pobj$method, direction = direction)
          
        } else {
          psfmi_boot <- psfmi_coxr(formula = pobj$formula_final, data = x, nimp=pobj$nimp, impvar = pobj$impvar,
                                   p.crit = p.crit, keep.predictors = pobj$keep.predictors,
                                   method = pobj$method, direction = direction)
        }
      }
      
      if(direction=="BW")
        boot_predictors_in <- ifelse(psfmi_boot$predictors_out[nrow(psfmi_boot$predictors_out), ], 0, 1)
      if(direction=="FW")
        boot_predictors_in <- psfmi_boot$predictors_in[nrow(psfmi_boot$predictors_in), ]
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
    bif_pat_sort <- data.frame(bif_pat %>%
                                 arrange(desc(n)) %>%
                                 select(1:ncol(bif_pat)))
    bif_pat_perc <- round((bif_pat_sort$n / nboot) * 100, 0)
    bif_pat_sort <- data.frame(bif_pat_sort, bif_pat_perc)
    
    if(!start_model) {
      colnames(bif_pat_sort) <- c(pobj$predictors_final, "freq", "bif_pat_perc")
    } else {
      colnames(bif_pat_sort) <- c(pobj$predictors_initial, "freq", "bif_pat_perc")
    }
    rownames(bif) <- paste("boot", 1:nboot)
    
    bif_total <- colSums(bif)
    bif_perc <- round((bif_total / nboot) * 100, 3)
    
    stabobj <- list(bif = bif, bif_total = bif_total, bif_perc = bif_perc,
                    model_stab = bif_pat_sort, call = call)
    return(stabobj)
  }
  
  if(inherits(pobj, "pmods") & boot_method == "cluster")
  {
    if(is_empty(pobj$clusvar))
      stop("\n", "No cluster variable defined, use single bootstrapping method")
    cluster_bootmi <- function(data, impvar, clusvar, nimp, nboot)
    {
      boot_imp <- list()
      for(i in 1:nimp){
        
        # Extract each imputed dataset
        imp_data <- subset(data, data[impvar] == i)
        nest_data <- imp_data %>% nest(data=c(-which(names(imp_data) == clusvar)))
        nest_data
        
        bs_data <- bootstraps(nest_data, times = nboot)
        bs_data
        
        df_boot <-lapply(bs_data$splits,
                         function(x) {
                           data_boot <- as_tibble(x) %>%
                             unnest(cols = c(clusvar))
                           n_row <- map_int(data_boot$data, ~ sapply(.x[1], NROW))
                           df_booti <- as_tibble(data.frame(clus_id=rep(data_boot[[clusvar]], n_row),
                                                            do.call("rbind", data_boot$data)))
                           names(df_booti)[1] <- clusvar
                           return(df_booti)
                         })
        boot_imp[[i]] <- df_boot
      }
      
      # binds first elements of lists
      dfs_bootmi <- transpose(boot_imp) %>% map(bind_rows)
      
      bootnr <- (1:nboot)
      nrow_boot <- map_int(dfs_bootmi, ~ sapply(.x[1], NROW))
      dfs_bootmi_tot <- mapply(function(x, y, z) {
        idboot <- rep(x, y)
        data.frame(idboot, z)},
        x=bootnr, y=nrow_boot, z=dfs_bootmi, SIMPLIFY = FALSE)
      return(dfs_bootmi_tot)
    }
    
    boot_clus_res <- cluster_bootmi(data=pobj$data, impvar = pobj$impvar,
                                    clusvar = pobj$clusvar, nimp = pobj$nimp, nboot = nboot)
    
    boot_pred_pat <- mapply(function(x, y) {
      message("\n", "Boot ", y)
      x <- as.data.frame(x)
      psfmi_boot <- psfmi_mm(data=x, nimp=pobj$nimp, impvar = pobj$impvar, random.eff = pobj$random.eff,
                             Outcome = pobj$Outcome, predictors = pobj$predictors, family = pobj$family,
                             p.crit = pobj$p.crit, cat.predictors = pobj$cat.predictors,
                             spline.predictors = pobj$spline.predictors, clusvar = pobj$clusvar,
                             int.predictors = pobj$int.predictors, keep.predictors = pobj$keep.predictors,
                             nknots = pobj$nknots, method = pobj$method, print.method = pobj$print.method)
      
      boot_predictors_in <- psfmi_boot$predictors_in[nrow(psfmi_boot$predictors_in), ]
      return(boot_predictors_in)
    }, x = boot_clus_res, y=boot_seq, SIMPLIFY = FALSE)
    
    
    bif <- data.frame(do.call("rbind", boot_pred_pat))
    
    if(!start_model)
      colnames(bif) <- pobj$predictors_final
    if(start_model)
      colnames(bif) <- pobj$predictors_initial
    
    names_temp <- colnames(bif)
    # Group selected models
    bif_pat <- bif %>%
      group_by_all() #%>%
    n <- count(bif_pat)$n
    bif_pat <- count(bif_pat)
    bif_pat <- data.frame(bif_pat, n)[, -length(colnames(bif_pat))]
    colnames(bif_pat) <- c(names_temp, "n")
    #names_temp <- names(bif_pat)
    
    bif_pat_sort <- arrange(bif_pat, desc(n))
    #bif_pat_sort <- data.frame(bif_pat %>%
    #                             arrange(desc(n)) %>%
    #                             select(1:ncol(bif_pat)))
    bif_pat_perc <- round((bif_pat_sort$n / nboot) * 100, 0)
    bif_pat_sort <- data.frame(bif_pat_sort, bif_pat_perc)
    
    bif_pat_perc <- round((bif_pat_sort$n / nboot) * 100, 0)
    bif_pat_sort <- data.frame(bif_pat_sort, bif_pat_perc)
    colnames(bif_pat_sort) <- c(names(pobj$predictors_in), "freq", "bif_pat_perc")
    
    rownames(bif) <- paste("boot", 1:nboot)
    
    bif_total <- colSums(bif)
    bif_perc <- round((bif_total / nboot) * 100, 3)
    
    stabobj <- list(bif = bif, bif_total = bif_total, bif_perc = bif_perc,
                    model_stab = bif_pat_sort, call = call)
    return(stabobj)
  }
}