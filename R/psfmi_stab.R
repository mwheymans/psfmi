#' Function to evaluate bootstrap predictor and model stability in multiply imputed datasets.
#'
#' \code{psfmi_stab} Stability analysis of predictors and prediction models selected with
#'  the \code{psfmi_lr}, \code{psfmi_coxr} or \code{psfmi_mm} functions of the \code{psfmi} package.
#'
#' @param pobj An object of class \code{smodsmi} (selected models in multiply imputed datasets), 
#'  produced by a previous call to \code{psfmi_lr}, \code{psfmi_coxr} or \code{psfmi_mm}.
#' @param boot_method A single string to define the bootstrap method. Use "single" after a call to 
#'  \code{psfmi_lr} and \code{psfmi_coxr} and "cluster" after a call to \code{psfmi_mm}.    
#' @param nboot A numerical scalar. Number of bootstrap samples to evaluate the stability. Default is 20.
#'
#' @details The function evaluates predictor selection frequency in stratified or cluster bootstrap samples.
#'  The stratification factor is the variable that separates the imputed datasets. It uses as input
#'  an object of class \code{smodsmi} as a result of a previous call to the \code{psfmi_lr}, 
#'  \code{psfmi_coxr} or \code{psfmi_mm} functions. In combination with the \code{psfmi_mm} function
#'  a cluster bootstrap method is used where bootstrapping is used on the level of the clusters only.  
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
#' @examples
#'  pool_lr <- psfmi_lr(data=lbpmilr, nimp=5, impvar="Impnr", Outcome="Chronic",
#'                    predictors=c("Gender", "Smoking",  "JobControl", "JobDemands",
#'                    "Age", "Radiation", "SocialSupport", "Function"),
#'                    cat.predictors = c("Carrying"), p.crit =0.157, method="D1")
#'  pool_lr$RR_Model
#'  pool_lr$multiparm
#'
#'  stab_res <- psfmi_stab(pool_lr, boot_method = "single", nboot=10) 
#'  stab_res$bif
#'  stab_res$bif_perc
#'  stab_res$model_stab
#'
#' @export
psfmi_stab <- function(pobj, boot_method=NULL, nboot=20)
{
  
  if(class(pobj)!="smodsmi") 
    stop("Object should be of type smodsmi", "\n")
  if(is.null(boot_method)) 
    stop("boot_method is not defined, choose single or cluster", "\n")
  if(pobj$p.crit==1) 
    stop("To determine Model Stability p.crit must be < 1", "\n")
  call <- match.call()
  nboot <- nboot
  data <- pobj$data
  boot_seq <- as.list(1:nboot)
  
  if (class(pobj)=="smodsmi" & boot_method == "single"){ 
    boot_data <- bootstraps(data, strata = pobj$impvar, times = nboot)
    boot_pred_pat <- mapply(function(x, y) {
      message("\n", "Boot ", y)
      x <- as.data.frame(x)
      if(pobj$model_type=="binomial") {
        psfmi_boot <- psfmi_lr(data=x, nimp=pobj$nimp, impvar = pobj$impvar,
                               Outcome = pobj$Outcome, predictors = pobj$predictors,
                               p.crit = pobj$p.crit, cat.predictors = pobj$cat.predictors,
                               spline.predictors = pobj$spline.predictors,
                               int.predictors = pobj$int.predictors, keep.predictors = pobj$keep.predictors,
                               knots = pobj$knots, method = pobj$method, print.method = pobj$print.method)
      }
      if(pobj$model_type=="survival") {
        psfmi_boot <- psfmi_coxr(data=x, nimp=pobj$nimp, impvar = pobj$impvar,
                                 time = pobj$time, status = pobj$status, predictors = pobj$predictors,
                                 p.crit = pobj$p.crit, cat.predictors = pobj$cat.predictors,
                                 spline.predictors = pobj$spline.predictors,
                                 int.predictors = pobj$int.predictors, keep.predictors = pobj$keep.predictors,
                                 knots = pobj$knots, method = pobj$method, print.method = pobj$print.method)
      }
      
      boot_predictors_in <- psfmi_boot$predictors_in[nrow(psfmi_boot$predictors_in), ]
      return(boot_predictors_in)
    }, x = boot_data$splits, y=boot_seq, SIMPLIFY = FALSE)
  }
  if(class(pobj)=="smodsmi" & boot_method == "cluster")
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
                             knots = pobj$knots, method = pobj$method, print.method = pobj$print.method)
      
      boot_predictors_in <- psfmi_boot$predictors_in[nrow(psfmi_boot$predictors_in), ]
      return(boot_predictors_in)
    }, x = boot_clus_res, y=boot_seq, SIMPLIFY = FALSE)
  }
  
  bif <- data.frame(do.call("rbind", boot_pred_pat))
  colnames(bif) <- names(pobj$predictors_in)
  
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
  colnames(bif_pat_sort) <- c(names(pobj$predictors_in), "freq", "bif_pat_perc")
  
  rownames(bif) <- paste("boot", 1:nboot)
  
  bif_total <- colSums(bif)
  bif_perc <- round((bif_total / nboot) * 100, 3)
  
  stabobj <- list(bif = bif, bif_total = bif_total, bif_perc = bif_perc,
                  model_stab = bif_pat_sort, call = call)
  return(stabobj)
}