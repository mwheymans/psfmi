#' Predictor selection function for Cox regression models
#'
#' \code{psfmi_coxr} Backward selection of predictors for Cox regression
#'  models using different selection methods.
#'
#' @param data Data frame or data matrix with stacked multiple imputed datasets.
#'   The original dataset that contains missing values must be excluded from the dataset.
#' @param nimp A numerical scalar. Number of imputed datasets. Default is 5.
#' @param impvar A character vector. Name of the variable that distinguishes the imputed datasets.
#' @param time Follow up time.
#' @param status The status variable, normally 0=censoring, 1=event.
#' @param predictors Character vector with the names of the predictor variables.
#'   At least one predictor variable has to be defined.
#' @param p.crit A numerical scalar. P-value selection criterium.
#' @param cat.predictors A single string or a vector of strings to define the categorical variables.
#'   Default is NULL categorical predictors.
#' @param int.predictors A single string or a vector of strings with the names of the variables that form
#'   an interaction pair, separated by a “:” symbol.
#' @param keep.predictors A single string or a vector of strings including the variables that are forced
#'   in the model during predictor selection. Categorical and interaction variables are allowed.
#' @param method A character vector to indicate the pooling method for p-values to pool the
#'   total model or used during predictor selection. This can be "D1", "D2", "MR" or "MPR".
#'   See details for more information.
#' @param print.method logical vector. If TRUE full matrix with p-values of all variables according to
#'   chosen method (under method) is shown. If FALSE (default) p-value for categorical variables according
#'   to method are shown and for continuous and dichotomous predictors Rubin’s Rules are used.
#'
#' @details print.method allows to choose between different p-value pooling methods. “D1” is pooling of
#'   the total covariance matrix, ”D2” is pooling of Chi-square values, “MR” is pooling Likelihood ratio
#'   statistics (method of Meng and Rubin) and “MPR” is pooling of median p-values (MPR rule).
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
#' @references Van Buuren S. (2012). Flexible Imputation of Missing Data. Chapman & Hall/CRC
#'   Interdisciplinary Statistics. Boca Raton.
#'
#' @examples
#'   # Pooling model (without backward selection) using D1
#'   psfmi_coxr(data=lbpmicox, nimp=5, impvar="Impnr", time="Time", status="Status",
#'   predictors=c("Duration", "Radiation", "Onset", "Function", "Age",
#'   "Previous", "Tampascale", "JobControl", "JobDemand", "Social"), p.crit=1,
#'   method="D1", cat.predictors=c("Expect_cat") )
#'
#'   # Predictor selection using p<0.05 and method D1
#'   psfmi_coxr(data=lbpmicox, nimp=5, impvar="Impnr", time="Time", status="Status",
#'   predictors=c("Duration", "Radiation", "Onset", "Function", "Age",
#'   "Previous", "Tampascale", "JobControl", "JobDemand", "Social"), p.crit=0.05,
#'   method="D1", cat.predictors=c("Expect_cat"))
#'
#'   # Predictor selection, force variable Tampa scale in model
#'   # using p<0.05 and method D1
#'   psfmi_coxr(data=lbpmicox, nimp=5, impvar="Impnr", time="Time", status="Status",
#'   predictors=c("Duration", "Previous",  "Radiation", "Onset",
#'   "Function", "Tampascale" ), p.crit=0.05, cat.predictors=c("Satisfaction",
#'   "Expect_cat"), int.predictors=c("Tampascale:Radiation",
#'   "Expect_cat:Tampascale"), keep.predictors = "Tampascale", method="D1" )
#'
#' @export
psfmi_coxr <-
  function(data, nimp=5, impvar=NULL, time, status, predictors=NULL,
  p.crit=1, cat.predictors=NULL, int.predictors=NULL, keep.predictors=NULL,
  method=NULL, print.method=FALSE)
{

    P <- predictors
    cat.P <- cat.predictors
    keep.P <- keep.predictors
    int.P <- int.predictors
    P.check <-c(P, cat.P)

    # Check data input
    if (!(is.matrix(data) | is.data.frame(data)))
      stop("Data should be a matrix or data frame")
    data <- as.data.frame(data)
    if ((nvar <- ncol(data)) < 2)
      stop("Data should contain at least two columns")
    if(is.null(impvar))
      stop("Imputation variable is not defined")
    if(is.null(method))
      stop("Define selection method: D1, D2, MR or MPR")
    if (order(unique(data[, impvar]))[1] == 0)
      stop("Original dataset should not be included")
    if(is.null(nimp))
      stop("Number of imputed datasets is not defined, use nimp!")
    if (nimp < 2) {
      stop("\n", "Number of imputed datasets must be > 1", "\n\n")
    }
    if (p.crit > 1) stop("\n", "P-value criterium > 1", "\n")
    if (!is.null(cat.P)) {
      if(any(cat.P%in%P)){
        cat.P.double <- cat.P[cat.P%in%P]
        cat(red("\n", "Categorical variable(s) -", cat.P.double,
          "- also defined as Predictor", "\n\n"))
        stop()
      }
    }

    if(any(duplicated(P))){
      cat(red("\n", "Predictor(s) - ", c(P[duplicated(P)]),
        " - defined more than once", "\n\n"))
      stop()
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
      stop("\n", "No predictors defined, cannot f model", "\n\n")

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

    coef.f <- se.f <- list()
    for (k in 1:length(P)) {
      cat("\n", "Step", k, "\n")
      cat("\n", "Variables included in model =", P, "\n")

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

    #### Rubin's Rules

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

    cat("\n", "Pooled model (Rubin's Rules)", "\n")
    print(pool.RR)

    if(any(is.infinite(unlist(pool.RR)))){
      cat(red("\n", "Check Pooled Model,
        some estimates could not be estimated", "\n"))
      stop()
    }

    ### D2
    if(method=="D2")
    {
    #### Get Chi/square values
    LChisq <- apply(chi.LR, 1 , as.list)

    mi.chiL <- lapply(1:length(LChisq),
      function(i, x, y) {
      x <- unlist(LChisq[[i]])
      y <- df.chi[[i]]
      miceadds::micombine.chisquare(x, y, display = F) })

    mi.chisq <- round(data.frame(do.call("rbind",
      mi.chiL))[, -c(3,4)], 4)

    rownames(mi.chisq) <- P
    if(print.method) {
      print(mi.chisq)
      cat("\n", "D2 Pooled p-values", "\n")
    }
    # Combine D2 with RR
    id.p.RR <- grep("factor", row.names(pool.RR))
    if(length(id.p.RR)==0){
      mi.chisq[, 2] <- pool.RR[, 3]
    }
    res.RR <- pool.RR[-id.p.RR, 3]
    mi.chisq[names(res.RR), 2] <- res.RR
    cat("\n", "Pooled p-values (D2 & RR)", "\n")
    print(mi.chisq)
    # Res.f is the dataframe with the CHISQ
    # pooled p-values for all variables
    }
    ### D1
    if(method=="D1")
    {
    est.D1 <- data.frame(D1_cox(data = data,
      impvar = impvar,
      nimp = nimp,
      fm = fm,
      names.var = names.var))
    rownames(est.D1) <- P
    if(print.method) {
      cat("\n", "D1 Pooled p-values", "\n")
      print(est.D1)
    }
    # Combine D1 with RR
    id.p.RR <- grep("factor", row.names(pool.RR))
    if(length(id.p.RR)==0){
        est.D1[, 2] <- pool.RR[, 3]
    }
    res.RR <- pool.RR[-id.p.RR, 3]
    est.D1[names(res.RR), 2] <- res.RR
    cat("\n", "Pooled p-values (D1 & RR)", "\n")
    print(est.D1)
    }
    ### Med P Rule
    if(method=="MPR")
    {
    med.pvalue <- round(data.frame(apply(chi.p,
       1, median)), 4)
    names(med.pvalue) <- "Median P-values"
    rownames(med.pvalue) <- P
     if(print.method) {
        cat("\n", "MPR Pooled p-values", "\n")
        print(med.pvalue)
     }
    # Combine Median p with RR
    id.p.RR <- grep("factor", row.names(pool.RR))
    if(length(id.p.RR)==0){
      med.pvalue[, 1] <- pool.RR[, 3]
    }
    res.RR <- pool.RR[-id.p.RR, 3]
    med.pvalue[names(res.RR), 1] <- res.RR
    cat("\n", "Pooled p-values (MPR & RR)", "\n")
    print(med.pvalue)
    }

    if(method=="D2"){
    p.pool <- data.frame(as.matrix(mi.chisq)[, 2])
    names(p.pool) <- "RR & D2 p-value"
    }
    if(method=="D1"){
    p.pool <- data.frame(as.matrix(est.D1)[, 2])
    names(p.pool) <- "RR & D1 p-value"
    }
    if(method=="MPR"){
    p.pool <- med.pvalue
    names(p.pool) <- "RR & MPR p-value"
    }

    if(!is.null(int.P)) {
      if(max(p.pool[, 1]) > p.crit) {
        del.coef.id.1 <-
          which(p.pool[, 1] == max(p.pool[, 1]))
        coef.excl.1 <- P[del.coef.id.1]
        if(p.crit < 1)
          cat("\n", "Variable with highest p-value at Step",
            k, "is", coef.excl.1, "\n")
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
      cat("\n", "Selection terminated. Variable(s) ",
          keep.P, "last variable(s) in model", "\n")
      (break)()
    }
   }
    del.coef.id <- which(p.pool[, 1] == max(p.pool[, 1]))
    if(length(del.coef.id)>1) {
      cat("\n", "Predictors with exact same P-value,
        first one excluded", "\n")
      del.coef.id <- del.coef.id[1]
    }
    coef.excl <- P[del.coef.id]
    if (p.pool[, 1][del.coef.id] > p.crit) {
      if(!is.null(keep.P)) cat("\n", "Predictor(s) to keep are",
                 keep.P.orig, "\n")
      cat("\n", "Variable excluded at Step",
          k, "is ", coef.excl, "\n")
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
        cat("\n", "Selection terminated. Variable(s) to keep -",
        keep.P.orig, " - in model", "\n")
      else cat("\n", "Selection terminated. No Variable(s) in model", "\n")
      (break)()
    }

    # Select variable with highest p-value
    del.coef.id <- which(p.pool[, 1] == max(p.pool[, 1]))
    if(length(del.coef.id)>1) {
      cat("\n", "Predictors with exact same P-value,
        first one excluded", "\n")
      del.coef.id <- del.coef.id[1]
    }
    coef.excl <- P[del.coef.id]
    if (p.pool[, 1][del.coef.id] > p.crit) {
      if(!is.null(keep.P)) cat("\n", "Predictor(s) to keep are",
          keep.P.orig, "\n")
    cat("\n", "Variable excluded at Step", k, "is -", coef.excl)
    cat("\n", "(Possibly, taking interaction terms or
      predictors to keep into account)", "\n")
    }
    P.drop <- grep(coef.excl, P.start, fixed=T)
    P <- P.start[-P.drop]
    if (method=="D1") names.var <- names.var.start[-P.drop]
    }

    if (p.pool[, 1][del.coef.id] > p.crit) {
      if (length(P) == 0) {
        cat("\n", "All variables excluded, computation terminated, ",
            "\n", "\n", "Model is empty", "\n")
        (break)()
      }
    Y <- c(paste("survival::Surv(", time, ",", status, ")~"))
    fm <- as.formula(paste(Y, paste(P, collapse = "+")))
    }

    if (p.pool[, 1][del.coef.id] < p.crit) {
      if(p.crit==1) cat("\n", "Final results of Pooled model
        with a p-value of", p.crit, "\n")
      if(!is.null(keep.P)){
        if(p.crit < 1) cat("\n", "Final result of Pooled model
          with a p-value of", p.crit,
              "and predictors to keep", keep.P.orig, "\n")
       else {
        cat("\n", "Final result of Pooled model
          with a p-value of", p.crit, "\n")
       }
      }
      break()
    }
  }
}