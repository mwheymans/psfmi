#' Pooling and Predictor selection function for multilevel
#' models in multiply imputed datasets
#'
#' \code{psfmi_mm} Pooling and backward selection for 2 level (generalized)
#' linear mixed models in multiply imputed datasets using different selection methods.
#'
#' @param data Data frame with stacked multiple imputed datasets.
#'   The original dataset that contains missing values must be excluded from the
#'   dataset. The imputed datasets must be distinguished by an imputation variable,
#'   specified under impvar, and starting by 1 and the clusters should be 
#'   distinguished by a cluster variable, specified under clusvar.
#' @param nimp A numerical scalar. Number of imputed datasets. Default is 5.
#' @param impvar A character vector. Name of the variable that distinguishes the
#'   imputed datasets.
#' @param clusvar A character vector. Name of the variable that distinguishes the
#'  clusters.
#' @param Outcome Character vector containing the name of the outcome variable.
#' @param predictors Character vector with the names of the predictor variables.
#'   At least one predictor variable has to be defined.
#' @param random.eff Character vector to specify the random effects as used by the 
#'   \code{lmer} and \code{glmer} functions of the \code{lme4} package.  
#' @param family Character vector to specify the type of model, "linear" is used to 
#'   call the \code{lmer} function and "binomial" is used to call the \code{glmer}
#'   function of the \code{lme4} package. See details for more information.
#' @param p.crit A numerical scalar. P-value selection criterium. A value of 1 
#'   provides the pooled model without selection.
#' @param cat.predictors A single string or a vector of strings to define the
#' categorical variables. Default is NULL categorical predictors.
#' @param spline.predictors A single string or a vector of strings to define the
#' (restricted cubic) spline variables. Default is NULL spline predictors. See details.
#' @param int.predictors A single string or a vector of strings with the names of the 
#' variables that form an interaction pair, separated by a “:” symbol.
#' @param keep.predictors A single string or a vector of strings including the variables that are forced
#'   in the model during predictor selection. Categorical and interaction variables are allowed.
#' @param knots A numerical vector that defines the number of knots for each spline predictor separately.
#' @param method A character vector to indicate the pooling method for p-values to pool the
#'   total model or used during predictor selection. This can be "D1", "D2", "D3" or "MPR".
#'   See details for more information.
#' @param print.method logical vector. If TRUE full matrix with p-values of all variables according to
#'   chosen method (under method) is shown. If FALSE (default) p-value for categorical variables according
#'   to method are shown and for continuous and dichotomous predictors Rubin’s Rules are used.
#'
#' @details The basic pooling procedure to derive pooled coefficients, standard errors, 95
#'  confidence intervals and p-values is Rubin's Rules (RR). Specific procedures are
#'  available to derive pooled p-values for categorical (> 2 categories) and spline variables.
#'  print.method allows to choose between the pooling methods: D1, D2 and D3 and MPR for pooling of 
#'  median p-values (MPR rule). The D1, D2 and D3 methods are called from the package \code{mitml}.
#'  For Logistic multilevel models (that are estimated using the \code{glmer} function), the D3 method
#'  is not yet available. Spline regression coefficients are defined by using the rcs function for 
#'  restricted cubic splines of the rms package. A minimum number of 3 knots as defined under knots is required.
#'
#'@return An object of class \code{smodsmi} (selected models in multiply imputed datasets) from 
#'  which the following objects can be extracted: imputed datasets as \code{data}, selected 
#'  pooled model as \code{RR_model}, pooled p-values according to pooling method as \code{multiparm}, 
#'  random effects as \code{random.eff}, predictors included at each selection step as \code{predictors_in}, 
#'  predictors excluded at each step as \code{predictors_out}, and \code{family}, \code{impvar}, \code{clusvar}, 
#'  \code{nimp}, \code{Outcome}, \code{method}, \code{p.crit}, \code{predictors}, \code{cat.predictors}, 
#'  \code{keep.predictors}, \code{int.predictors}, \code{spline.predictors}, \code{knots}, \code{print.method}, 
#'  \code{model_type} and \code{call} .
#'
#' @references Eekhout I, van de Wiel MA, Heymans MW. Methods for significance testing of categorical
#'   covariates in logistic regression models after multiple imputation: power and applicability
#'   analysis. BMC Med Res Methodol. 2017;17(1):129.
#' @references Enders CK (2010). Applied missing data analysis. New York: The Guilford Press.
#' @references Meng X-L, Rubin DB. Performing likelihood ratio tests with multiply-imputed data sets.
#'   Biometrika.1992;79:103-11.
#' @references van de Wiel MA, Berkhof J, van Wieringen WN. Testing the prediction error difference between
#'   2 predictors. Biostatistics. 2009;10:550-60.
#' @references mitml package https://cran.r-project.org/web/packages/mitml/index.html
#' @references Van Buuren S. (2018). Flexible Imputation of Missing Data. 2nd Edition. Chapman & Hall/CRC
#'   Interdisciplinary Statistics. Boca Raton.
#' @references http://missingdatasolutions.rbind.io/
#'
#' @examples
#'   pool_mm <- psfmi_mm(data=ipdna_md, nimp=5, impvar=".imp", family="linear",
#'   predictors=c("gender", "afib", "sbp"), clusvar = "centre",
#'   random.eff="( 1 | centre)", Outcome="dbp", cat.predictors = "bmi_cat",
#'   p.crit=0.15, method="D1", print.method = FALSE)
#'   pool_mm$RR_Model
#'   pool_mm$multiparm
#'
#' @export
psfmi_mm <- function(data, nimp=5, impvar=NULL, clusvar = NULL, Outcome, predictors=NULL, 
                     random.eff=NULL, family="linear", p.crit=1, cat.predictors=NULL, 
                     spline.predictors=NULL, int.predictors=NULL, keep.predictors=NULL, 
                     knots=NULL, method="RR", print.method=FALSE)
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
if(family == "binomial" & !all(data[Outcome]==1 | data[Outcome]==0))
  stop("Outcome should be a 0 - 1 variable")
if ((nvar <- ncol(data)) < 2)
  stop("Data should contain at least two columns")
if(is.null(impvar))
  stop("Imputation variable is not defined")
if(is.null(clusvar))
  stop("Cluster variable is not defined")
if(is.null(random.eff))
  stop("Random effects not specified")
if(is.null(method))
  stop("Define selection method: D1, D2, D3 or MPR")
if (order(unique(data[, impvar]))[1] == 0)
  stop("Original dataset should not be included")
if(is.null(nimp))
  stop("Number of imputed datasets is not defined, use nimp!")
if (nimp < 2) {
  stop("\n", "Number of imputed datasets must be > 1", "\n\n")
}
if (p.crit > 1)
  stop("\n", "P-value criterium > 1", "\n")
if (any(knots < 3))
  stop("\n", "Number of knots must be > 2", "\n")
if (length(knots) != length(s.P))
  stop("\n", "Number of knots not specified for every spline variable", "\n")
if (!is.null(cat.P)) {
  if(any(cat.P %in% P)){
    cat.P.double <- cat.P[cat.P %in% P]
    stop("\n", "Categorical variable(s) -", cat.P.double,
         "- also defined as Predictor", "\n\n")
  }
}
if (!is.null(s.P)){
  if(any(s.P %in% P)){
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
# First predictors, second cetegorical
# predictors and last interactions
P <- c(P, cat.P, s.P, int.P)
if (is.null(P))
  stop("\n", "No predictors defined, cannot fit model", "\n\n")
# Define predictors from model 
# order of interaction term changes
if (!is.null(int.P)) {
  Y <- c(paste(Outcome, paste("~")))
  fm <- as.formula(paste(Y, paste(c(P, random.eff), collapse = "+")))
  f <- lmer(fm, data = data[data[impvar] == 1, ], REML = FALSE) 
  if(family=="binomial") f <- glmer(fm, data = data[data[impvar] == 1, ],
               family = binomial)

P <- names(fixef(f))[-1]
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
  } else {
    for(i in 1:length(cat.P)) {
      P <- gsub(cat.P[i],
                replacement=paste0("factor(", cat.P[i], ")"), P)
      if(!is.null(keep.P)){
        keep.P <- gsub(cat.P[i],
                       replacement=paste0("factor(", cat.P[i], ")"), keep.P)
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

keep.P.orig <- keep.P

if(any(!keep.P %in% P))
  stop("\n", "Variables to keep not defined as Predictor", "\n\n")

# Start  loop for backward selection over imputed datasets

coef.f <- se.f <- RR.model <- multiparm <- coef.excl_step <- step.nr <- P_in_step <- list()

for (k in 1:length(P)) {
  P_in_step[[k]] <- P
  
  if(method=="D3" & family=="binomial"){
    message("\n", "The D3 method is not available yet for Logistic multilevel models", "\n",
            "Use D1, D2 or MPR instead", "\n")
  }
  
    chi.LR <- data.frame(matrix(0, length(P), nimp))
    chi.p <- data.frame(matrix(0, length(P), nimp))
    
    Y <- c(paste(Outcome, paste("~")))
    fm <- as.formula(paste(Y, paste(c(P, random.eff), collapse = "+")))
    
    # Extract df of freedom for pooling Chi-Square values
    df.chi <- as.list(car::Anova(lmer(fm, REML = FALSE,
               data=data[data[impvar] == 1, ]))$'Df')
    if(family=="binomial") df.chi <- as.list(car::Anova(glmer(fm,
             data=data[data[impvar] == 1, ], family=binomial))$'Df')
    
    # Start loop for Rubin's Rules
    for (i in 1:nimp) {
      
      f <- lmer(fm, data = data[data[impvar] == i, ], REML = FALSE)
      if(family=="binomial") f <- glmer(fm, data = data[data[impvar] == i, ],
                   family = binomial)

      coef.f[[i]] <- summary(f)[[10]][, 1]
      se.f[[i]] <- summary(f)[[10]][, 2]
      chi.LR[, i] <- car::Anova(f)$`Chisq`
      chi.p[, i] <- car::Anova(f)$`Pr(>Chisq)`
    }
    coef.f.qhat <- do.call("rbind", coef.f)
    
    # Rubin's Rules
    RR <- norm::mi.inference(coef.f, se.f, 0.95)
    pool.RR <- do.call("cbind", RR)[, -c(3, 7, 8)]
    
    if(family=="binomial"){
      OR <- exp(pool.RR[, 1])
      L.OR <- exp(pool.RR[, 4])
      U.OR <- exp(pool.RR[, 5])
      pool.RR <- round(cbind(pool.RR, OR, L.OR, U.OR), 5)
    } else{
      pool.RR <- round(pool.RR, 5)
    }
    
    RR.model[[k]] <- pool.RR
    names(RR.model)[k] <- paste("Step", k)
    
    #if(any(is.infinite(pool.RR))){
    #  stop("\n", "Check Pooled Model, some parameters
    #      could not be estimated", "\n")
    #}
    p.pool <- data.frame(pool.RR[-1, 3])
    if(method=="RR") multiparm <- NULL
    
    if(family=="linear" | (family=="binomial" & method!="D3")) {
    
      # D1
      if(method=="D1" | method=="D2" | method=="D3") {
        est.DX <- est.DX_orig <- data.frame(psfmi_mm_multiparm(data = data, nimp = nimp, impvar = impvar,
                Outcome = Outcome, P = P, family = family, random.eff = random.eff,  
                method = method, print.method = print.method))
        rownames(est.DX) <- rownames(est.DX_orig)  <- P  
        
        # Combine DX with RR
        id.p.RR.f <- grep("factor", row.names(pool.RR))
        id.p.RR.spl <- grep("rcs", row.names(pool.RR))
        res.RR <- pool.RR[-c(1, id.p.RR.f,id.p.RR.spl), 3]
        est.DX[names(res.RR), 1] <- res.RR
        names(est.DX) <- c( paste(method, "and RR p-values"), paste(method, "F-statistic"))
       
        if(print.method) {
          est.DX <- est.DX_orig
          names(est.DX) <- c( paste(method, "p-values"), paste(method, "F-statistic"))
        }
        
        multiparm[[k]] <- est.DX
        names(multiparm)[k] <- paste("Step", k)
      }
    # MPR
    if(method=="MPR") {
        med.pvalue <- med.pvalue_orig <- round(data.frame(apply(chi.p, 1, median)), 5)
        rownames(med.pvalue) <- rownames(med.pvalue_orig) <- P
        
        # Combine Median p with RR
        id.p.RR.f <- grep("factor", row.names(pool.RR))
        id.p.RR.spl <- grep("rcs", row.names(pool.RR))
        res.RR <- pool.RR[-c(1, id.p.RR.f, id.p.RR.spl), 3]
        med.pvalue[names(res.RR), 1] <- res.RR
        names(med.pvalue) <- ("MPR & RR P-values")
        
        if(print.method) {
          med.pvalue <- med.pvalue_orig
          names(med.pvalue) <- ("MPR P-values")
        }
        
        multiparm[[k]] <- med.pvalue
        names(multiparm)[k] <- paste("Step", k)
      }
    
      if(method=="D1" | method=="D2" | method=="D3") p.pool <- est.DX
      if(method=="MPR") p.pool <- med.pvalue

}
    
    if(!is.null(int.P)) {
      if(max(p.pool[, 1]) > p.crit) {
        del.coef.id.1 <- which(p.pool[, 1] == max(p.pool[, 1]))
        coef.excl.1 <- P[del.coef.id.1]
      }
    }
    
    if(!any(grepl(":", P))) {
      P.start <- P
      
      if (!is.null(keep.P)) {
        id.var.exclude <- lapply(keep.P, function(x) {
          grep(x, P, fixed=T)
        })
        
        P <- P[-unlist(id.var.exclude)]
        
        p.pool <- data.frame(p.pool[-c(unlist(id.var.exclude)), ])
        names(p.pool) <- NULL
        
        if (nrow(p.pool) == 0) {
          message("\n", " Selection correctly terminated. Variable(s) ",
                  paste(keep.P, collapse = " ") , " last variable(s) in model", "\n")
          (break)()
        }
      }
      
      del.coef.id <- which(p.pool[, 1] == max(p.pool[, 1]))
      if(length(del.coef.id)>1) {
        #cat("\n", "Predictors with exact same P-value,
        #  first one chosen", "\n")
        del.coef.id <- del.coef.id[1]
      }
      coef.excl <- P[del.coef.id]
      if (p.pool[, 1][del.coef.id] > p.crit) {
        message("Variable excluded at Step ",
                k, " is - ", coef.excl)
      }
      P.drop <- grep(coef.excl, P.start, fixed=T)
      P <- P.start[-P.drop]
      
    }
    
    # Identify variables part of interaction term
    # Define when interaction term will be deleted
    if(any(grepl(":", P))) {
      P.start <- P
      
      i.int.var <- lapply(P[grep(":", P)], function(x) {
        unlist(strsplit(x, split=":"))
      })
      i.int.var <- unique(unlist(i.int.var))
      
      # Evaluate if there are more interaction terms in model and
      # if variables that are part of interaction term
      # are also part of other interaction terms and if so,
      # exclude variables before selection
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
      P <- P[-id.var.exclude]
      
      # Determine if interaction term has to
      # be excluded before predictor selection
      p.pool <- data.frame(p.pool[-c(id.var.exclude), ],
                           row.names = P)
      names(p.pool) <- NULL
      
      if (!is.null(keep.P)) {
        id.var.keep <- lapply(keep.P, function(x) {
          grep(x, P, fixed=T)
        })
        
        P <- P[-unlist(id.var.keep)]
        
        # Determine if keep.Predictors must be
        # excluded before predictor selection (id.var.exclude)
        p.pool <- data.frame(p.pool[-c(unlist(id.var.keep)), ],
                             row.names = P)
        names(p.pool) <- NULL
      }
      
      if (nrow(p.pool) == 0) {
        if(!is.null(keep.P)) message("\n", "Selection correctly terminated.
          Variable(s) to keep - ", paste(keep.P.orig, collapse=" "),
                                     " - last variable(s) in model", "\n")
        else message("\n", "Model is empty after last step", "\n")
        (break)()
      }
      
      # Select variable with highest p-value
      del.coef.id <- which(p.pool[, 1] == max(p.pool[, 1]))
      if(length(del.coef.id)>1) {
        #cat("\n", "Predictors with exact same P-value,
        #  first one chosen", "\n")
        del.coef.id <- del.coef.id[1]
      }
      coef.excl <- P[del.coef.id]
      
      if (p.pool[, 1][del.coef.id] > p.crit) {
        message("Variable excluded at Step ", k,
                " is - ", coef.excl)
      }
      P.drop <- grep(coef.excl, P.start, fixed=T)
      P <- P.start[-P.drop]
      
    }
    
    if (p.pool[, 1][del.coef.id] > p.crit) {
      if (length(P) == 0) {
        coef.excl_step[[k]] <- coef.excl
        message("\n", "Selection correctly terminated, ",
                "\n", "Model is empty after last step", "\n")
        (break)()
      }
      Y <- c(paste(Outcome, paste("~")))
      fm <- as.formula(paste(Y, paste(c(P, random.eff), collapse = "+")))
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
  P_select <- P_in_step[[1]]
}
else {
  if(is_empty(coef.excl_step)){
    coef.excl_step <- as.null(coef.excl_step)
    P_select <- data.frame(matrix(rep(1, length(P_in_step[[1]])), 1, 
                                  length(P_in_step[[1]]), byrow=TRUE))
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
}

pobj <- list(data = data, RR_Model = RR.model, multiparm = multiparm, random.eff = random.eff,
                  predictors_in = P_select, predictors_out = coef.excl_step, family = family,
                  impvar = impvar, clusvar = clusvar, nimp = nimp, Outcome = Outcome, method = method, p.crit = p.crit,
                  predictors = predictors, cat.predictors = cat.predictors, 
                  keep.predictors = keep.predictors, int.predictors = int.predictors, model_type = family,
                  spline.predictors = spline.predictors, knots = knots, print.method = print.method, call = call)
class(pobj) <- "smodsmi"
return(pobj)
}