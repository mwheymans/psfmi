#' Pooling and predictor selection function for Cox regression
#' models in multiply imputed datasets
#'
#' \code{psfmi_coxr} Pooling and backward selection for Cox regression
#' models in multiply imputed datasets using different selection methods.
#'
#' @param data Data frame with stacked multiple imputed
#' datasets. The original dataset that contains missing values must be excluded from the dataset. The imputed
#'   datasets must be distinguished by an imputation variable, specified under impvar, and starting by 1.
#' @param nimp A numerical scalar. Number of imputed datasets. Default is 5.
#' @param impvar A character vector. Name of the variable that distinguishes the imputed datasets.
#' @param time Follow up time.
#' @param status The status variable, normally 0=censoring, 1=event.
#' @param predictors Character vector with the names of the predictor variables.
#'   At least one predictor variable has to be defined.
#' @param p.crit A numerical scalar. P-value selection criterium.
#' @param cat.predictors A single string or a vector of strings to define the categorical variables.
#'   Default is NULL categorical predictors.
#' @param spline.predictors A single string or a vector of strings to define the
#' (restricted cubic) spline variables. Default is NULL spline predictors.
#' @param int.predictors A single string or a vector of strings with the names of the variables that form
#'   an interaction pair, separated by a “:” symbol.
#' @param keep.predictors A single string or a vector of strings including the variables that are forced
#'   in the model during predictor selection. Categorical and interaction variables are allowed. See details.
#' @param knots A numerical vector that defines the number of knots for each spline predictor separately.
#' @param method A character vector to indicate the pooling method for p-values to pool the
#'   total model or used during predictor selection. This can be "D1", "D2", or "MPR".
#'   See details for more information.
#' @param print.method logical vector. If TRUE full matrix with p-values of all variables according to
#'   chosen method (under method) is shown. If FALSE (default) p-value for categorical variables according
#'   to method are shown and for continuous and dichotomous predictors Rubin’s Rules are used.
#'
#' @details The basic pooling procedure to derive pooled coefficients, standard errors, 95
#'  confidence intervals and p-values is Rubin's Rules (RR). Specific procedures are
#'  available to derive pooled p-values for categorical (> 2 categories) and spline variables.
#'  print.method allows to choose between these pooling methods that are:
#'  “D1” is pooling of the total covariance matrix, ”D2” is pooling of Chi-square values,
#'  and “MPR” is pooling of median p-values (MPR rule). Spline regression coefficients are defined
#'  by using the rcs function for restricted cubic splines of the rms package of Frank Harrell.
#'  A minimum number of 3 knots as defined under knots is needed.
#'
#'@return An object of class \code{smodsmi} (selected models in multiply imputed datasets) from 
#'  which the following objects can be extracted: imputed datasets as \code{data}, selected 
#'  pooled model as \code{RR_model}, pooled p-values according to pooling method as \code{multiparm}, 
#'  predictors included at each selection step as \code{predictors_in}, predictors excluded at each step 
#'  as \code{predictors_out}, and \code{impvar}, \code{nimp}, \code{time}, \code{status}, \code{method}, 
#'  \code{p.crit}, \code{predictors}, \code{cat.predictors}, \code{keep.predictors}, \code{int.predictors}, 
#'  \code{spline.predictors}, \code{knots}, \code{print.method}, \code{call}, \code{model_type} and 
#'  \code{predictors_final} for names of predictors in final selection step, \code{fit.formula} is the 
#'  regression formula of start model and \code{predictors_initial} for names of predictors in start
#'  model.
#'    
#' @references Eekhout I, van de Wiel MA, Heymans MW. Methods for significance testing of categorical
#'   covariates in logistic regression models after multiple imputation: power and applicability
#'   analysis. BMC Med Res Methodol. 2017;17(1):129.
#' @references Enders CK (2010). Applied missing data analysis. New York: The Guilford Press.
#' @references van de Wiel MA, Berkhof J, van Wieringen WN. Testing the prediction error difference between
#'   2 predictors. Biostatistics. 2009;10:550-60.
#' @references Marshall A, Altman DG, Holder RL, Royston P. Combining estimates of interest in prognostic
#'   modelling studies after multiple imputation: current practice and guidelines. BMC Med Res Methodol.
#'   2009;9:57.
#' @references Van Buuren S. (2018). Flexible Imputation of Missing Data. 2nd Edition. Chapman & Hall/CRC
#'   Interdisciplinary Statistics. Boca Raton.
#' @references http://missingdatasolutions.rbind.io/
#' 
#' @examples
#'   pool_coxr <- psfmi_coxr(data=lbpmicox, nimp=5, impvar="Impnr", time="Time", 
#'   status="Status", predictors=c("Duration", "Radiation", "Onset"), p.crit=1, 
#'   method="D1", cat.predictors=c("Expect_cat"))
#'   pool_coxr$RR_Model
#'   pool_coxr$multiparm
#'   
#'   pool_coxr <- psfmi_coxr(data=lbpmicox, nimp=5, impvar="Impnr", time="Time", 
#'   status="Status", predictors=c("Previous",  "Radiation", "Onset",
#'   "Function", "Tampascale" ), p.crit=0.05, cat.predictors=c("Expect_cat"), 
#'   int.predictors=c("Tampascale:Radiation",
#'   "Expect_cat:Tampascale"), keep.predictors = "Tampascale", method="D2")
#'   pool_coxr$RR_Model
#'   pool_coxr$multiparm
#'   pool_coxr$predictors_in
#'   
#' @export
psfmi_coxr <-
  function(data, nimp=5, impvar=NULL, time, status, predictors=NULL,
   p.crit=1, cat.predictors=NULL, spline.predictors=NULL, int.predictors=NULL,
   keep.predictors=NULL, knots=NULL, method="RR", print.method=FALSE)
  {
    call <- match.call()
    
    P <- predictors
    cat.P <- cat.predictors
    keep.P <- keep.predictors
    int.P <- int.predictors
    s.P <- spline.predictors
    P.check <-c(P, cat.P, s.P)
    
    # Check data input
    if (!(is.data.frame(data)))
      stop("Data should be a data frame")
    data <- data.frame(as_tibble(data))
    data <- mutate_if(data, is.factor, ~ as.numeric(as.character(.x)))
    if(!all(data[status]==1 | data[status]==0))
      stop("Outcome should be a 0 - 1 variable")
    if ((nvar <- ncol(data)) < 2)
      stop("Data should contain at least two columns")
    if(is.null(impvar))
      stop("Imputation variable is not defined")
    if(is.null(method)) method="RR"
    if(all(!is.null(cat.predictors) | !is.null(spline.predictors)) & method=="RR")
      stop("Categorical or spline variables in model, define selection method: D1, D2, D3 or MPR")
    if (order(unique(data[, impvar]))[1] == 0)
      stop("Original dataset should not be included")
    if(is.null(nimp))
      stop("Number of imputed datasets is not defined, use nimp!")
    if (nimp < 2) {
      stop("\n", "Number of imputed datasets must be > 1", "\n\n")
    }
    if (p.crit > 1) stop("\n", "P-value criterium > 1", "\n")
    if (any(knots<3))
      stop("\n", "Number of knots must be > 2", "\n")
    if (length(knots) != length(s.P))
      stop("\n", "Number of knots not specified for every spline variable", "\n")
    if (!is.null(cat.P)) {
      if(any(cat.P%in%P)){
        cat.P.double <- cat.P[cat.P%in%P]
        stop("\n", "Categorical variable(s) -", cat.P.double,
                "- also defined as Predictor", "\n\n")
      }
    }
    if (!is.null(s.P)) {
      if(any(s.P%in%P)){
        s.P.double <- s.P[s.P%in%P]
        stop("\n", "Spline variable(s) -", s.P.double,
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
    
    # First predictors, second categorical
    # predictors and last interactions!
    P <- c(P, cat.P, s.P, int.P)
    if (is.null(P))
      stop("\n", "No predictors defined, cannot fit model", "\n\n")
    
    # Define predictors from model for D1 method because order
    # of variables of interaction term changes
    if (!is.null(int.P)) {
      Y <- c(paste("survival::Surv(", time, ",", status, ")~"))
      fm <- as.formula(paste(Y, paste(P, collapse = "+")))
      f <- survival::coxph(fm, data = data[data[impvar] == 1, ])
      P <- names(coef(f))
    }
    if (!is.null(keep.P)) {
      for(i in 1:length(keep.P)){
        if(grepl(":", keep.P[i])) {
          keep.P.spl <- unlist(strsplit(keep.P[i], split=":"))
          if(length(P[Reduce("&", lapply(keep.P.spl, grepl, P))])==0)
            stop("Interaction term in keep.predictors not defined
              as int.predictors, incorrect")
          keep.P[i] <- P[Reduce("&", lapply(keep.P.spl, grepl, P))]
        }
      }
    }
    
    if (!is.null(cat.P)) {
      if(length(cat.P)==1){
        P <- gsub(cat.P,
           replacement=paste0("factor(", cat.P, ")"), P)
        if(!is.null(keep.P)){
          keep.P <- gsub(cat.P,
           replacement=paste0("factor(", cat.P, ")"), keep.P)
        }
      }
      else {
        for(i in 1:length(cat.P)){
          P <- gsub(cat.P[i], replacement=paste0("factor(",
           cat.P[i], ")"), P)
          if(!is.null(keep.P)){
            keep.P <- gsub(cat.P[i], replacement=paste0("factor(",
              cat.P[i], ")"), keep.P)
          }
        }
      }
    }
    if (!is.null(s.P)) {
      if(length(s.P)==1){
        P <- gsub(s.P,
                  replacement=paste0("rcs(", s.P, ",", knots, ")"), P)
        if(!is.null(keep.P)){
          keep.P <- gsub(s.P,
                         replacement=paste0("rcs(", s.P, ",", knots, ")"), keep.P)
        }
      } else {
        for(i in 1:length(s.P)) {
          P <- gsub(s.P[i],
                    replacement=paste0("rcs(", s.P[i], ",", knots[i], ")"), P)
          if(!is.null(keep.P)){
            keep.P <- gsub(s.P[i],
                           replacement=paste0("rcs(", s.P[i], ",", knots[i], ")"), keep.P)
          }
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
    
    if(method=="D1") {
      names.var <- list()
      if(is.null(int.P)){
        for(i in 1:length(P)){
          Y <- c(paste("survival::Surv(", time, ",", status, ")~"))
          fm <- as.formula(paste(Y, paste(P[i], collapse = "+")))
          f <- survival::coxph(fm, data = data[data[impvar] == 1, ])
          if(any(is.na(f$coefficients))){
            print(summary(f))
            stop("Startmodel too complex. Contains singularities")
          }
          if(grepl("factor", P[[i]])) {
            names.var[[i]] <- names(coef(f))
          }
          else names.var[[i]] <- names(coef(f))
        }
      }
      if(!is.null(int.P)){
        for(i in 1:length(P)){
          P[i] <- gsub(":", "*", P[i])
          # Create names.var for D1 method and for complicated
          # combinations as interaction between categorical variables
          Y <- c(paste("survival::Surv(", time, ",", status, ")~"))
          fm <- as.formula(paste(Y, paste(P[i], collapse = "+")))
          f <- survival::coxph(fm, data = data[data[impvar] == 1, ])
          if(any(is.na(f$coefficients))){
            print(summary(f))
            stop("Startmodel too complex. Contains singularities")
          }
          if(grepl("[*]", P[[i]])){
            names.var[[i]] <- names(coef(f))
            names.var[[i]] <- names.var[[i]][grepl(":",
                                                   names.var[[i]])]
          }
          else names.var[[i]] <- names(coef(f))
          P[i] <- gsub("[*]", ":", P[i])
        }
      }
    }
    keep.P.orig <- keep.P
    
    if(any(!keep.P %in% P))
      stop("\n", "Variables to keep not defined as Predictor", "\n\n")
    
    coef.f <- se.f <- RR.model <- multiparm <- coef.excl_step <- step.nr <- P_in_step <- list()
    for (k in 1:length(P)) {
      P_in_step[[k]] <- P
      chi.LR <- data.frame(matrix(0,
                                  length(P), nimp))
      chi.p <- data.frame(matrix(0,
                                 length(P), nimp))
      
      Y <- c(paste("survival::Surv(", time, ",", status, ")~"))
      fm <- as.formula(paste(Y, paste(P, collapse = "+")))
      
      # Extract df of freedom for pooling Chi-Square values
      if (length(P) == 1) {
        df.chi <- as.list(car::Anova(coxph(fm,
             data=data[data[impvar] == 1, ]))$'Df'[2])
      }
      if (length(P) > 1) {
        df.chi <- as.list(car::Anova(coxph(fm,
             data=data[data[impvar] == 1, ]))$'Df')
      }
      
      for (i in 1:nimp) {
        f <- survival::coxph(fm, data = data[data[impvar] == i, ])
        coef.f[[i]] <- summary(f)[[7]][, 1]
        names(coef.f[[i]]) <- names(coef(f))
        se.f[[i]] <- summary(f)[[7]][, 3]
        names(se.f[[i]]) <- names(coef(f))
        if (length(P) == 1) {
          chi.p[i] <- car::Anova(f)$Pr[2]
          chi.LR[i] <- car::Anova(f)$'Chisq'[2]
        }
        if (length(P) > 1) {
          chi.p[, i] <- car::Anova(f)$`Pr(>Chisq)`
          chi.LR[, i] <- car::Anova(f)$'LR Chisq'
        }
      }
      
      coef.f.qhat <- do.call("rbind", coef.f)
      
      # Rubin's Rules
      
      RR <- norm::mi.inference(coef.f, se.f, 0.95)
      pool.RR <- do.call("cbind", RR)[, -c(3, 7, 8)]
      pool.RR <- round(pool.RR, 4)
      
      if(length(RR$est)==1){
        pool.RR <- as.data.frame(t(do.call("cbind",
          RR)[, -c(3, 7, 8)]))
        row.names(pool.RR) <- P
        pool.RR <- round(pool.RR, 4)
      }
      
      HR <- exp(pool.RR[, 1])
      L.HR <- exp(pool.RR[, 4])
      U.HR <- exp(pool.RR[, 5])
      pool.RR <- round(cbind(pool.RR, HR,
                             L.HR, U.HR), 4)
      
      RR.model[[k]] <- pool.RR
      names(RR.model)[k] <- paste("Step", k)
      
      if(any(is.infinite(unlist(pool.RR)))){
        stop("\n", "Check Pooled Model,
        some parameters could not be estimated", "\n")
      }
      if(method=="RR"){
        p.pool <- data.frame(pool.RR[, 3])
        multiparm <- NULL
      }
      # D2
      if(method=="D2")
      {
        # Get Chi/square values
        LChisq <- apply(chi.LR, 1 , as.list)
        
        mi.chiL <- lapply(1:length(LChisq),
           function(i, x, y) {
           x <- unlist(LChisq[[i]])
           y <- df.chi[[i]]
           miceadds::micombine.chisquare(x, y, display = F) })
        
        mi.chisq <- mi.chisq_orig <- round(data.frame(do.call("rbind",
                    mi.chiL))[, -c(3,4)], 4)
        row.names(mi.chisq) <- row.names(mi.chisq_orig) <- P
        
        # Combine D2 with RR
        id.p.RR.f <- grep("factor", row.names(pool.RR))
        id.p.RR.spl <- grep("rcs", row.names(pool.RR))
        res.RR <- pool.RR[-c(id.p.RR.f, id.p.RR.spl), 3]
        mi.chisq[names(res.RR), 2] <- res.RR
        names(mi.chisq) <- c("D2", "p-value D2 & RR")
        
        if(print.method) {
          mi.chisq <- mi.chisq_orig
          names(mi.chisq) <- c("D2", "p-value D2")
        }
        
        multiparm[[k]] <- mi.chisq
        names(multiparm)[k] <- paste("Step", k)
        # Res.f is the dataframe with the CHISQ
        # pooled p-values for all variables
      }
      # D1
      if(method=="D1")
      {
        est.D1 <- est.D1_orig <- data.frame(D1_cox(data = data,
          impvar = impvar,
          nimp = nimp,
          fm = fm,
          names.var = names.var))
        row.names(est.D1) <- row.names(est.D1_orig) <- P
        
        # Combine D1 with RR
        id.p.RR.f <- grep("factor", row.names(pool.RR))
        id.p.RR.spl <- grep("rcs", row.names(pool.RR))
        res.RR <- pool.RR[-c(id.p.RR.f, id.p.RR.spl), 3]
        est.D1[names(res.RR), 2] <- res.RR
        names(est.D1) <- c("Chi_sq", "p-value D1 & RR")
        
        if(print.method) {
          est.D1 <- est.D1_orig
          names(est.D1) <- c("Chi_sq", "p-value D1")
        }
        
        multiparm[[k]] <- est.D1
        names(multiparm)[k] <- paste("Step", k)
      }
      # Med P Rule
      if(method=="MPR")
      {
        med.pvalue <- med.pvalue_orig <- round(data.frame(apply(chi.p,
                                                                1, median)), 4)
        row.names(med.pvalue) <- row.names(med.pvalue_orig) <- P
        names(med.pvalue) <- "P-value MPR & RR"
        
        # Combine Median p with RR
        id.p.RR.f <- grep("factor", row.names(pool.RR))
        id.p.RR.spl <- grep("rcs", row.names(pool.RR))
        res.RR <- pool.RR[-c(id.p.RR.f, id.p.RR.spl), 3]
        med.pvalue[names(res.RR), 1] <- res.RR
        
        if(print.method) {
          med.pvalue <- med.pvalue_orig
          names(med.pvalue) <- "P-value MPR"
        }
        
        multiparm[[k]] <- med.pvalue
        names(multiparm)[k] <- paste("Step", k)
      }
      
      if(method=="D2"){
        p.pool <- data.frame(as.matrix(mi.chisq)[, 2])
        names(p.pool) <- "p-value RR & D2"
      }
      if(method=="D1"){
        p.pool <- data.frame(as.matrix(est.D1)[, 2])
        names(p.pool) <- "p-value RR & D1"
      }
      if(method=="MPR"){
        p.pool <- med.pvalue
        names(p.pool) <- "p-value RR & MPR"
      }
      
      if(!is.null(int.P)) {
        if(max(p.pool[, 1]) > p.crit) {
          del.coef.id.1 <-
            which(p.pool[, 1] == max(p.pool[, 1]))
          coef.excl.1 <- P[del.coef.id.1]
        }
      }
      
      if(!any(grepl(":", P))) {
        P.start <- P
        if (method=="D1")
          names.var.start <- names.var
        
        if (!is.null(keep.P)) {
          id.var.exclude <- lapply(keep.P,
              function(x) { grep(x, P, fixed=T)
              })
          
          if (method=="D1")
            names.var <- names.var[-unlist(id.var.exclude)]
          P <- P[-unlist(id.var.exclude)]
          
          p.pool <- data.frame(p.pool[-c(unlist(id.var.exclude)), ])
          names(p.pool) <- NULL
          
          if (nrow(p.pool) == 0) {
            message("\n", "Selection correctly terminated. Variable(s) ",
              paste(keep.P, collapse = " "), " last variable(s) in model", "\n")
            (break)()
          }
        }
        del.coef.id <- which(p.pool[, 1] == max(p.pool[, 1]))
        if(length(del.coef.id)>1) {
          #  cat("\n", "Predictors with exact same P-value,
          #first one excluded", "\n")
          del.coef.id <- del.coef.id[1]
        }
        coef.excl <- P[del.coef.id]
        if (p.pool[, 1][del.coef.id] > p.crit) {
          message("Variable excluded at Step ", k,
                  " is - ", coef.excl)
        }
        P.drop <- grep(coef.excl, P.start, fixed=T)
        P <- P.start[-P.drop]
        if (method=="D1") names.var <- names.var.start[-P.drop]
      }
      
      # Identify which variables are part of interaction term
      # Define when interaction term will be deleted
      # Interaction with categorical variable
      if(any(grepl(":", P))) {
        P.start <- P
        if (method=="D1") names.var.start <- names.var
        
        i.int.var <- lapply(P[grep(":", P)], function(x) {
          unlist(strsplit(x, split=":"))
        })
        i.int.var <- unique(unlist(i.int.var))
        
        # Evaluate if there are more interaction terms in model and
        # if variables that are part of interaction term
        # are also part of other interaction terms and if so,
        # exclude all variables before selection
        int.var.terms <- grep(":", P)
        if(length(int.var.terms)>1){
          names.int.var.terms <- P[int.var.terms ]
          i.int.var <- lapply(names.int.var.terms[grep(":",
             names.int.var.terms)], function(x) {
             unlist(strsplit(x, split=":"))
             })
          i.int.var <- unique(unlist(i.int.var))
        }
        
        id.int.var.exl <- lapply(i.int.var, function(x) {
          grep(x, P, fixed = T)
        })
        id.int.var.exl <- unique(unlist(id.int.var.exl))
        
        names.int.id <- grep(":", P[id.int.var.exl])
        id.var.exclude <- id.int.var.exl[-names.int.id]
        
        P.var.int.keep <- P[id.var.exclude]
        if (method=="D1") names.var <- names.var[-id.var.exclude]
        P <- P[-id.var.exclude]
        
        # Determine if interaction term has to
        # be excluded before predictor selection (id.var.exclude)
        p.pool <- data.frame(p.pool[-c(id.var.exclude), ], row.names = P)
        names(p.pool) <- NULL
        
        if (!is.null(keep.P)) {
          id.var.keep <- lapply(keep.P, function(x) {
            grep(x, P, fixed=T)
          })
          
          if (method=="D1") names.var <- names.var[-unlist(id.var.keep)]
          P <- P[-unlist(id.var.keep)]
          
          # Determine if keep.predictors must be
          # excluded before predictor selection (id.var.exclude2)
          p.pool <- data.frame(p.pool[-c(unlist(id.var.keep)), ],
                               row.names = P)
          names(p.pool) <- NULL
        }
        
        if (nrow(p.pool) == 0) {
          if(!is.null(keep.P))
            message("\n", "Selection correctly terminated. Variable(s) to keep - ",
             paste(keep.P.orig, collapse=" "), " - last variable(s) in model", "\n")
          else message("\n", "Model is empty after last step", "\n")
          (break)()
        }
        
        # Select variable with highest p-value
        del.coef.id <- which(p.pool[, 1] == max(p.pool[, 1]))
        if(length(del.coef.id)>1) {
          #cat("\n", "Predictors with exact same P-value,
          #first one excluded", "\n")
          del.coef.id <- del.coef.id[1]
        }
        coef.excl <- P[del.coef.id]
        if (p.pool[, 1][del.coef.id] > p.crit) {
          message("Variable excluded at Step ", k, " is - ", coef.excl)
        }
        P.drop <- grep(coef.excl, P.start, fixed=T)
        P <- P.start[-P.drop]
        if (method=="D1") names.var <- names.var.start[-P.drop]
      }
      
      if (p.pool[, 1][del.coef.id] > p.crit) {
        if (length(P) == 0) {
          coef.excl_step[[k]] <- coef.excl
          message("\n", "Selection correctly terminated, ",
                  "\n", "Model is empty after last step", "\n")
          (break)()
        }
        Y <- c(paste("survival::Surv(", time, ",", status, ")~"))
        fm <- as.formula(paste(Y, paste(P, collapse = "+")))
      }
      
      if (p.pool[, 1][del.coef.id] < p.crit) {
        if(p.crit==1) message("\n", "Pooled model correctly estimated
            using a p-value of ", p.crit, "\n")
        if(!is.null(keep.P)){
          if(p.crit < 1) message("\n", "Pooled model correctly estimated
          using a p-value of ", p.crit,
          " and predictors to keep ", paste(keep.P.orig, collapse=" "), "\n")
        }
        break()
      }
      coef.excl_step[[k]] <- coef.excl
      # End k loop
    }
    if(p.crit==1) {
      coef.excl_step <- as.null(coef.excl_step)
      P_select <- data.frame(matrix(1, 1, length(P_in_step[[1]]), byrow=TRUE))
      colnames(P_select) <- P_in_step[[1]]
      predictors_final <- P_in_step[[1]]
    }
    else {
      if(is_empty(coef.excl_step)){
        coef.excl_step <- as.null(coef.excl_step)
        P_select <- data.frame(matrix(1, 1, length(P_in_step[[1]]), byrow=TRUE))
        rownames(P_select) <- "Step 1"
        colnames(P_select) <- P_in_step[[1]]
      }
      else{  
        coef.excl_step <- data.frame(do.call("rbind", coef.excl_step))
        names(coef.excl_step) <- "Excluded"
        row.names(coef.excl_step) <- paste("Step", 1:nrow(coef.excl_step))
        
        outOrder_step <- P_in_step[[1]]
        P_select <- data.frame(do.call("rbind", lapply(P_in_step, function(x) {
          outOrder_step %in% x
        })))
        names(P_select) <- P_in_step[[1]]
        P_select[P_select==TRUE] <- 1
        row.names(P_select) <- paste("Step", 1:nrow(P_select))
        if(length(P_in_step[[1]]) == length(coef.excl_step[, 1])) {
          r_null <- rep(0, length(P_in_step[[1]]))
          names(r_null) <- P_in_step[[1]]
          P_select <- bind_rows(P_select, r_null)
          row.names(P_select) <- paste("Step", 1:nrow(P_select))
        }
      }
      predictors_final <- colnames(P_select)[which(P_select[nrow(P_select), ]==1)]
    }
    fit.formula <- as.formula(paste(Y, paste(P_in_step[[1]], collapse = "+")))
    pobj <- list(data = data, RR_Model = RR.model, multiparm = multiparm,
                      predictors_in = P_select, predictors_out = coef.excl_step, time = time,
                      status = status, impvar = impvar, nimp = nimp, method = method, p.crit = p.crit,
                      predictors = predictors, cat.predictors = cat.predictors, call = call,
                      keep.predictors = keep.predictors, int.predictors = int.predictors, model_type = "survival",
                      spline.predictors = spline.predictors, knots = knots, print.method = print.method,
                      fit.formula = fit.formula, predictors_final = predictors_final,
                      predictors_initial = P_in_step[[1]])
    class(pobj) <- "smodsmi"
    return(pobj)
}
