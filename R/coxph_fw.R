#' Predictor selection function for forward selection of Cox regression models.
#'
#' \code{coxph_bw} Forward selection of Cox regression models in single dataset 
#'  using as selection method the partial likelihood-ratio statistic.
#'
#' @param data A data frame. 
#' @param formula A formula object to specify the model as normally used by coxph.
#'   See under "Details" and "Examples" how these can be specified.
#' @param time Survival time.
#' @param status The status variable, normally 0=censoring, 1=event.
#' @param predictors Character vector with the names of the predictor variables.
#'   At least one predictor variable has to be defined. Give predictors unique names
#'   and do not use predictor name combinations with numbers as, age2, gnder10, etc.
#' @param p.crit A numerical scalar. P-value selection criterium. A value of 1
#'   provides the pooled model without selection.
#' @param cat.predictors A single string or a vector of strings to define the
#' categorical variables. Default is NULL categorical predictors.
#' @param spline.predictors A single string or a vector of strings to define the
#' (restricted cubic) spline variables. Default is NULL spline predictors. See details.
#' @param int.predictors A single string or a vector of strings with the names of the variables that form
#'   an interaction pair, separated by a “:” symbol.
#' @param keep.predictors A single string or a vector of strings including the variables that are forced
#'   in the model during predictor selection. All type of variables are allowed.
#' @param nknots A numerical vector that defines the number of knots for each spline predictor separately.
#'
#' @details 
#'  A typical formula object has the form \code{Surv(time, status) ~ terms}. Categorical variables has to
#'  be defined as \code{Surv(time, status) ~ factor(variable)}, restricted cubic spline variables as
#'  \code{Surv(time, status) ~ rcs(variable, 3)}. Interaction terms can be defined as
#'  \code{Surv(time, status) ~ variable1*variable2} or \code{Surv(time, status) ~ variable1 + variable2 + 
#'  variable1:variable2}. All variables in the terms part have to be separated by a "+".
#'
#'@return An object of class \code{smods} (single models) from
#'  which the following objects can be extracted: original dataset as \code{data}, final selected
#'  model as \code{RR_model_final}, model at each selection step \code{RR_model},
#'  p-values at final step \code{multiparm_final}, and at each step as \code{multiparm}, 
#'  formula object at final step as \code{formula_final}, 
#'  and at each step as \code{formula_step} and for start model as \code{formula_initial}, 
#'  predictors included at each selection step as \code{predictors_in}, predictors excluded
#'  at each step as \code{predictors_out}, and \code{time}, \code{status}, \code{p.crit}, \code{call},
#'  \code{model_type}, \code{predictors_final} for names of predictors in final selection step and 
#'  \code{predictors_initial} for names of predictors in start model and \code{keep.predictors} for
#'  variables that are forced in the model during selection.
#'
#'@references http://missingdatasolutions.rbind.io/
#'
#'@examples
#' lbpmicox1 <- subset(psfmi::lbpmicox, Impnr==1) # extract first imputed dataset
#' res_single <- coxph_bw(data=lbpmicox1, p.crit = 0.05, formula=Surv(Time, Status) ~
#'                            Previous +  Radiation + Onset + Age + Tampascale + 
#'                            Pain + JobControl + factor(Satisfaction), 
#'                            spline.predictors = "Function",
#'                            nknots = 3)
#'          
#' res_single$RR_model_final
#' res_single$multiparm_final
#' 
#' @author Martijn Heymans, 2021
#' 
#' @export
coxph_fw <- function(data,
                     formula = NULL,
                     status = NULL,
                     time = NULL,
                     predictors=NULL,
                     p.crit=1,
                     cat.predictors=NULL,
                     spline.predictors=NULL,
                     int.predictors=NULL,
                     keep.predictors=NULL,
                     nknots=NULL)
{
  
  call <- match.call()
  
  if(is_empty(formula)) {
    if(!all(data[status]==1 | data[status]==0))
      stop("Status should be a 0 - 1 variable")
    P <-
      predictors
    cat.P <-
      cat.predictors
    int.P <-
      gsub(":", "*", int.predictors)
    s.P <-
      spline.predictors
  } else{
    form <-
      terms(formula)
    form_vars <-
      attr(form, "term.labels")
    if(is_empty(form_vars))
      stop("\n", "No predictors defined, model is empty")
    time <-
      as.character(attr(form, "variables")[[2]][[2]])
    status <-
      as.character(attr(form, "variables")[[2]][[3]])
    int.P <-
      form_vars[grepl(paste(c("[*]", ":"), collapse = "|"), form_vars)]
    int.P_temp <-
      unique(unlist(str_split(int.P, paste(c("[*]", ":"), collapse = "|"))))
    form_vars <-
      form_vars[!grepl(paste(c("[*]", ":"), collapse = "|"), form_vars)]
    form_vars <-
      unique(c(form_vars, int.P_temp))
    cat.P <-
      form_vars[grepl("factor", form_vars)]
    form_vars <-
      form_vars[!grepl("factor", form_vars)]
    s.P <-
      form_vars[grepl("rcs", form_vars)]
    nknots <-
      c(readr::parse_number(s.P))
    form_vars <-
      form_vars[!grepl("rcs", form_vars)]
    int.P <-
      gsub(":", "*", clean_P(int.P))
    cat.P <- clean_P(cat.P)
    s.P <- clean_P(s.P)
    P <- form_vars
  }
  keep.P <-
    gsub(":", "*", keep.predictors)
  keep.P <-
    sapply(as.list(keep.P), clean_P)
  
  P.check <-
    c(P, cat.P, s.P)
  # Check data input
  if(!(is.data.frame(data)))
    stop("Data should be a data frame")
  data <- data.frame(as_tibble(data))
  data <- mutate_if(data, is.factor, ~ as.numeric(as.character(.x)))
  if(!all(data[status]==1 | data[status]==0))
    stop("Status should be a 0 - 1 variable")
  if ((nvar <- ncol(data)) < 2)
    stop("Data should contain at least two columns")
  if (p.crit > 1)
    stop("\n", "P-value criterium > 1", "\n")
  if (any(nknots<3))
    stop("\n", "Number of knots must be > 2", "\n")
  if (length(nknots) != length(s.P))
    stop("\n", "Number of knots not specified for every spline variable", "\n")
  if (!is.null(cat.P)) {
    if(any(cat.P%in%P)){
      cat.P.double <- cat.P[cat.P%in%P]
      stop("\n", "Categorical variable(s) -", cat.P.double,
           "- also defined as Predictor", "\n\n")
    }
  }
  if (!is_empty(s.P)){
    if(any(s.P%in%P)){
      s.P.double <- s.P[s.P%in%P]
      stop("\n", "Do not include Spline variable(s) -", s.P.double,
           "- in predictors", "\n\n")
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
  if(!is_empty(int.P)) {
    int.P.check <- lapply(int.P[grep("[*]", int.P)],
                          function(x) { unlist(strsplit(x, split="[*]")) })
    int.P.check <- unique(unlist(int.P.check))
    if(any(!int.P.check %in% P.check))
      stop("\n", "Not all interaction terms defined as
        Predictor or Categorical Predictor", "\n\n")
  }
  # First predictors, second categorical
  # predictors and last interactions
  P <-
    c(P, cat.P, s.P, int.P)
  if (is.null(P))
    stop("\n", "No predictors to select, model is empty", "\n\n")
  
  if (!is_empty(keep.P)) {
    for(i in 1:length(keep.P)){
      if(grepl("[*]", keep.P[i])) {
        keep.P.spl <- unlist(strsplit(keep.P[i], split="[*]"))
        if(length(P[Reduce("&", lapply(keep.P.spl, grepl, P))])==0)
          stop("Interaction term in keep.predictors not defined
            as int.predictors, incorrect")
        keep.P[i] <- P[Reduce("&", lapply(keep.P.spl, grepl, P))]
      }
    }
  }
  
  if (!is_empty(cat.P)) {
    if(length(cat.P)==1){
      P <-
        gsub(cat.P,
             replacement=paste0("factor(", cat.P, ")"), P)
      if(!is.null(keep.P)){
        keep.P <-
          gsub(cat.P,
               replacement=paste0("factor(", cat.P, ")"), keep.P)
      }
    } else {
      for(i in 1:length(cat.P)) {
        P <-
          gsub(cat.P[i],
               replacement=paste0("factor(", cat.P[i], ")"), P)
        if(!is.null(keep.P)){
          keep.P <-
            gsub(cat.P[i],
                 replacement=paste0("factor(", cat.P[i], ")"), keep.P)
        }
      }
    }
  }
  if (!is_empty(s.P)) {
    if(length(s.P)==1){
      P <-
        gsub(s.P,
             replacement=paste0("rcs(", s.P, ",", nknots, ")"), P)
      if(!is.null(keep.P)){
        keep.P <-
          gsub(s.P,
               replacement=paste0("rcs(", s.P, ",", nknots, ")"), keep.P)
      }
    } else {
      for(i in 1:length(s.P)) {
        P <- gsub(s.P[i],
                  replacement=paste0("rcs(", s.P[i], ",", nknots[i], ")"), P)
        if(!is.null(keep.P)){
          keep.P <-
            gsub(s.P[i],
                 replacement=paste0("rcs(", s.P[i], ",", nknots[i], ")"), keep.P)
        }
      }
    }
  }
  levels.cat.P <- lapply(cat.P, function(x) {
    nr.levels.cat.P <- length(table(data[, x]))
    if (nr.levels.cat.P < 3) {
      stop("\n", "Categorical variable(s) only 2 levels,
        do not define as categorical", "\n\n")
    }
  })
  
  
  if(any(!keep.P %in% P))
    stop("\n", "Variables to keep not defined as Predictor", "\n\n")
  
  P_orig <-
    P
  P_orig_temp <-
    clean_P(P)
  keep.P <-
    clean_P(keep.P)
  
  if(!is_empty(keep.P))
    if(length(P_orig)==1)
      if(P_orig_temp == keep.P)
        stop("\n", "No need to define keep.predictors. Exclude keep.predictors and set p.crit = 1","\n")
  
  ############################## BW selection

  P_each_step <- fm_step <- fm_total <- RR_model_total <- chi_test_total <-
    RR_model_select <- multiparm_step <- multiparm_end <- list()
  
  P_select <- 0
  fm_step <-  as.list(rep(0, length(P)))
  
  if(!is_empty(keep.P)){
    P_temp <- clean_P(P)
    keep.P <-
      sapply(as.list(keep.P), clean_P)
    if(any(grepl("[*]", keep.P))){
      keep.P <- c(unique(unlist(str_split(keep.P[grep("[*]",
                                                      keep.P)], "[*]"))), keep.P)
    }
    keep.P_temp <- P[which(P_temp %in% keep.P)]
    P <- P[-which(P_temp %in% keep.P)]
    if(is_empty(P)) P <- keep.P_temp
    keep.P <- keep.P_temp
  }
  
  # Start J loop, to build up models,
  # variable by variable
  for(j in 1:length(P))
  {
    
    P_Chi <- chi_test_step <- RR.model <-
      fm_step <- multiparm <- as.list(rep(0, length(P)))
    
    # Loop k
    for (k in 1:length(P)) {
      
      # set regression formula fm
      Y <-
        c(paste0("Surv(", time, ",", status, ")~"))
      fm <- terms(as.formula(paste0(Y, paste0(c(P[k], keep.P), collapse = "+"))))
      
      if(P_select!=0){
        fm <- terms(update.formula(fm, paste0("~. +",
                                              paste0(paste0(P_each_step, collapse = "+")))))
      }
      
      if(P_select==0){
        cov.nam0 <- "1"
        cov.nam0_int <- cov.nam0_keep <- NULL
        # Test interaction terms against separate main effects
        if(!is_empty(keep.P))
          cov.nam0_keep <- keep.P
        if(grepl("[*]", P[[k]]))
          cov.nam0_int <- unlist(str_split(P[[k]], "[*]"))
        cov.nam0_temp <- unique(c(cov.nam0_keep, cov.nam0_int))
        if(!is_empty(cov.nam0_temp))
          cov.nam0 <- cov.nam0_temp
        f0 <- as.formula(paste(Y, paste(c(cov.nam0), collapse = "+")))
      }
      if(P_select!=0){
        f0_orig <- terms(update.formula(f0, paste0("~. +",
                                                   paste0(paste0(P_select, collapse = "+")))))
        cov.nam0_int <- NULL
        if(grepl("[*]", P[[k]]))
          cov.nam0_int <- unlist(str_split(P[[k]], "[*]"))
        f0 <- terms(update.formula(f0, paste0("~. +",
                                              paste0(paste0(c(P_select, cov.nam0_int), collapse = "+")))))
      }
      
      fit1 <- coxph(formula(fm), data = data)
      fit0 <- coxph(formula(f0), data = data)
      
      # Model estimates
      out.res <-
        summary(fit1)$coefficients
      lower.EXP <-
        exp(exp(out.res[, 1]) - (1.96*out.res[, 2]))
      upper.EXP <-
        exp(exp(out.res[, 1]) + (1.96*out.res[, 2]))
      model.res <-
        data.frame(cbind(out.res, lower.EXP, upper.EXP))
      names(model.res) <- c("Estimate", "HR", "Std Error", "Z value",
                            "P-value", "low EXP(HR)", "High EXP(HR)")
      
      RR.model[[k]] <- model.res
      names(RR.model)[k] <- paste("Step", j)
      if(P_select==0) names(RR.model)[k] <- paste("Step", 1)
      
      ll0 <-
        logLik(fit0)
      ll1 <-
        logLik(fit1)
      
      LL <-
        -2*(as.numeric(ll0) - as.numeric(ll1))
      df0 <-
        attr(ll0, "df")
      df1 <-
        attr(ll1, "df")
      diff_df <-
        df1 - df0
      pvalue <-
        pchisq(LL, df = diff_df, lower.tail = FALSE)
      
      pool.multiparm <- data.frame(matrix(c(pvalue, LL), length(P[k]), 2))
      row.names(pool.multiparm) <- P[k]
      names(pool.multiparm) <- c("p-values", "LR-statistic")
      pool.multiparm
      
      P_Chi[[k]] <- pool.multiparm
      
      # Set f0 to original, before testing interactions
      if(P_select==0){
        cov.nam0 <- "1"
        if(!is_empty(keep.P)) cov.nam0 <- keep.P
        f0 <- as.formula(paste(Y, paste(c(cov.nam0), collapse = "+")))
      }
      if(P_select!=0)
        f0 <- f0_orig
      # Extract regression formula's
      fm_step[[k]] <- paste(Y, paste(attr(fm, "term.labels"), collapse = " + "))
      names(fm_step)[k] <- paste("Test - ", P[k])
      
    }
    # End k loop
    ##############################################################
    
    fm_total[[j]] <- fm_step
    RR_model_total[[j]] <- RR.model
    
    # p.pool for multiparameer pooling
    p.pool <- data.frame(do.call("rbind", P_Chi))
    rownames_temp <- row.names(p.pool)
    p.pool <- data.frame(p.pool)#[, 1])
    row.names(p.pool) <- clean_P(rownames_temp)
    names(p.pool) <- c("p-value", "LR-statistic") #paste("p-value", "LR-statistic")
    
    if(P_select==0) names(fm_total)[[j]] <- paste("Step 0")
    else names(fm_total)[[j]] <- paste("Step", j-1)
    
    multiparm_end[[j]] <- p.pool
    # Extract variable with lowest P
    P_in <- which(p.pool[, 1] == min(p.pool[, 1]))
    if(length(P_in) > 1) {
      P_in <- P_in[1]
    }
    P_select <- P[P_in]
    
    # If selected predictor is interaction term
    # exclude separate variables from P list
    P_in_temp <- NULL
    P_select_temp <- row.names(p.pool)[P_in]
    if(grepl("[*]", P_select)){
      P_int_split <- unlist(str_split(P_select_temp, "[*]"))
      P_in_temp <- which(row.names(p.pool) %in% P_int_split)
    }
    
    if (p.pool[, 1][P_in] > p.crit) {
      message("\n", "Selection correctly terminated, ",
              "\n", "No new variables entered the model", "\n")
      P_each_step <- c(P_each_step[-j])
      if(is_empty(P_each_step)){
        fm_step <- as.formula(paste(Y, 1))
        if(!is_empty(keep.P)){
          fm_step <- as.formula(paste(Y, paste(keep.P, collapse = "+")))
          P_each_step <- c(keep.P)
        }
        #fit <- list()
        fit <- coxph(fm_step, data = data)
        
        # Model estimates
        model.res <- summary(fit)$coefficients
        RR_model_select <- list(model.res)
        names(RR_model_select)[[1]] <- paste("Step", 0, " - no variables entered - ")
        
        multiparm <- list(p.pool)
        names(multiparm)[[1]] <- paste("Step", 0, " - no variables entered - ")
      }
      (break)()
    }
    
    RR_model_select[[j]] <- RR.model[[P_in]]
    names(RR_model_select)[[j]] <- paste("Step", j, "- entered -", P_select)
    
    # Variables included in each step
    P_each_step[[j]] <- P_select
    
    if (p.pool[, 1][P_in] < p.crit) {
      message("Entered at Step ", j,
              " is - ", P_select)
    }
    
    row.names(p.pool) <- P
    P <- c(P[-c(P_in, P_in_temp)])
    multiparm_step[[j]] <- p.pool
    names(multiparm_step)[[j]] <- paste("Step", j-1, "- selected -", P_select)
    
    # P = 0, means all variables are included during FW selection
    if(is_empty(P)){
      message("\n", "Selection correctly terminated, ",
              "\n", "all variables added to the model", "\n")
      P_each_step <- c(P_each_step)#, keep.P)
      break()
    }
    # End J loop
  }
  
  # Extract selected models
  outOrder_step <- P_orig
  if(!is_empty(P_each_step)){
    P_select <- data.frame(do.call("rbind", lapply(P_each_step, function(x) {
      x <- str_replace(x, "[*]", ":")
      x <- unique(c(x, keep.P))
      outOrder_step %in% x
    })))
    names(P_select) <- P_orig
    row.names(P_select) <- paste("Step", 1:length(P_each_step))
    if(!nrow(P_select)==1) {
      P_select <- apply(P_select, 2, function(x) ifelse(x, 1, 0))
      P_select_final <- ifelse(colSums(P_select)>0, 1, 0)
      P_select <- rbind(P_select, P_select_final)
      row.names(P_select)[nrow(P_select)] <- "Included"
    } else {
      P_select_final <- P_select
      P_select <- rbind(P_select, P_select_final)
      P_select <- apply(P_select, 2, function(x) ifelse(x, 1, 0))
      row.names(P_select) <- c("Step 1", "Included")
    }
  } else {
    P_select <- matrix(rep(0, length(P_orig)), 1, length(P_orig))
    dimnames(P_select) <- list("Included", P_orig)
  }
  
  multiparm_out <- NULL
  if(length(c(P_select)==0)==1) { P_excluded <- P_orig
  } else {
    P_excluded <- as_tibble(names(P_select[nrow(P_select), ][P_select[nrow(P_select), ] ==0] ))
  }
  if(is_empty(P_excluded)){
    P_excluded <- NULL
  } else {
    names(P_excluded) <- "Excluded"
  }
  predictors_final <- names(P_select[nrow(P_select), ][P_select[nrow(P_select), ] ==1])
  
  if(is_empty(P_each_step)){
    RR_model <- RR_model_final <- RR_model_select
    multiparm <- multiparm
    fm_total <- fm_total
    names(RR_model_final) <- "Final model"
    multiparm_final <- multiparm
    names(multiparm_final) <- names(multiparm) <- "Step 0 - no variables entered"
    fm_step_final <- fm_total
    multiparm_out <- multiparm_end
    names(multiparm_out) <- "Predictors removed"
  }
  if(!is_empty(P_each_step)){
    RR_model <- RR_model_select
    multiparm <- multiparm_step
    fm_total <- fm_total
    if(is_empty(P))
    {
      RR_model_final <- RR_model[j]
      names(RR_model_final) <- "Final model"
      multiparm_final <- car::Anova(fit1, test.statistic="LR")
      fm_step_final <- fm_total[j]
    }
    if(!is_empty(P)){
      RR_model_final <- RR_model[j-1]
      multiparm_final <- multiparm[j-1]
      fm_step_final <- fm_total[j-1]
      if(j==1 & !is_empty(keep.P)) {
        RR_model_final <- RR_model
        fm_step_final <- fm_total
        multiparm <- car::Anova(fit, test.statistic="LR")
        multiparm_final <- car::Anova(fit, test.statistic="LR")
      }
      names(RR_model_final) <- "Final model"
      multiparm_out <- multiparm_end[j]
      names(multiparm_out) <- "Predictors removed"
    }
  }
  
  Y_initial <-
    c(paste0("Surv(", time, ",", status, ")~"))
  formula_initial <-
    as.formula(paste(Y_initial, paste(P_orig, collapse = "+")))
  
  pobjfw <-
    list(data = data, RR_model_final = RR_model_final, RR_model = RR_model,
         multiparm = multiparm, multiparm_final = multiparm_final,
         multiparm_out = multiparm_out, formula_step = fm_total,
         formula_final = fm_step_final, formula_initial = formula_initial,
         predictors_in = P_select, predictors_out = P_excluded, status = status,
         time = time, p.crit = p.crit, call = call, model_type = "survival",
         predictors_final = predictors_final, predictors_initial = P_orig,
         keep.predictors = keep.P)
  
  class(pobjfw) <- "smods"
  return(pobjfw)
}
