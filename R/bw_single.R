#' Predictor selection function for backward selection of
#' Logistic regression models.
#'
#' \code{bw_single} Backward selection of Logistic regression
#' prediction models using as selection method the likelihood-ratio Chi-square value.
#'
#' @param data A data frame. 
#' @param formula A formula object to specify the model as normally used by glm.
#'   See under "Details" and "Examples" how these can be specified.
#' @param Outcome Character vector containing the name of the outcome variable.
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
#'  A typical formula object has the form \code{Outcome ~ terms}. Categorical variables has to
#'  be defined as \code{Outcome ~ factor(variable)}, restricted cubic spline variables as
#'  \code{Outcome ~ rcs(variable, 3)}. Interaction terms can be defined as
#'  \code{Outcome ~ variable1*variable2} or \code{Outcome ~ variable1 + variable2 + variable1:variable2}.
#'  All variables in the terms part have to be separated by a "+".
#'
#'@return An object of class \code{smods} (single models) from
#'  which the following objects can be extracted: original dataset as \code{data}, final selected
#'  model as \code{RR_model_final}, model at each selection step \code{RR_model_setp},
#'  p-values at final step according to selection method as \code{multiparm_final}, and
#'  at each step as \code{multiparm_step}, formula object at final step as \code{formula_final}, 
#'  and at each step as \code{formula_step} and for start model as \code{formula_initial}, 
#'  predictors included at each selection step as \code{predictors_in}, predictors excluded
#'  at each step as \code{predictors_out}, and \code{Outcome}, \code{anova_test}, \code{p.crit}, \code{call},
#'  \code{model_type}, \code{predictors_final} for names of predictors in final selection step and 
#'  \code{predictors_initial} for names of predictors in start model.
#'
#'@references http://missingdatasolutions.rbind.io/
#'
#'@examples
#' res_single <- bw_single(data=lbpmilr, p.crit = 0.05, Outcome="Chronic",
#'          predictors=c("Tampascale", "Smoking"),
#'          cat.predictors = c("Satisfaction"))
#'          
#' res_single$RR_model_final
#'
#' @seealso \code{\link{psfmi_perform}}
#' 
#' @author Martijn Heymans, 2020
#' 
#' @export
bw_single <- function(data,
         formula = NULL,
         Outcome = NULL,
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
  if(is_empty(Outcome))
    stop("Outcome variable not defined")
  P <-
    predictors
  cat.P <-
    cat.predictors
  int.P <-
    int.predictors
  s.P <-
    spline.predictors
} else{
  form <-
    terms(formula)
  form_vars <-
    attr(form, "term.labels")
  if(is_empty(form_vars))
    stop("\n", "No predictors defined, model is empty")
  Outcome <-
    as.character(attr(form, "variables")[[2]])
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
if (!(is.data.frame(data)))
  stop("Data should be a data frame")
data <- data.frame(as_tibble(data))
data <- mutate_if(data, is.factor, ~ as.numeric(as.character(.x)))
if(!all(data[Outcome]==1 | data[Outcome]==0))
  stop("Outcome should be a 0 - 1 variable")
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

RR.model <- df.chi <- P_rm_step <- fm_step <- imp.dt <- multiparm <- list()

  # Loop k, to pool models in multiply imputed datasets
  for (k in 1:(length(P)+1)) {

    Y <-
      c(paste(Outcome, paste("~")))
    fm <-
      as.formula(paste(Y, paste(P, collapse = "+")))

      #fit <- list()
        fit <-
          glm(fm, data = data, family = binomial)
        chi_test <-
          car::Anova(fit)

      # Rubin's Rules
      out.res <-
        summary(fit)$coefficients
      OR <-
        exp(out.res[, 1])
      lower.EXP <-
        exp(exp(out.res[, 1]) - (1.96*out.res[, 2]))
      upper.EXP <-
        exp(exp(out.res[, 1]) + (1.96*out.res[, 2]))
      model.res <-
        data.frame(cbind(out.res, OR, lower.EXP, upper.EXP))
      names(model.res) <- c("Estimate", "Std Error", "Z value",
                            "P-value", "OR", "low EXP(OR)", "High EXP(OR)")
      RR.model[[k]] <-
        model.res
      names(RR.model)[[k]] <-
        paste("Step", k)#

      p.pool <-
        data.frame(matrix(0, length(P), 2))
      p.pool[, 1] <-
        chi_test$`Pr(>Chisq)`
      p.pool[, 2] <-
        chi_test$`LR Chisq`

      #P <-
      #  clean_P(P)
      row.names(p.pool) <-
        P
      names(p.pool) <-
        c("P-value", "LR Chisq")

      # Extract regression formula's
      fm_step[[k]] <-
        formula(fm)
      names(fm_step)[[k]] <-
        paste("Step", k)
      # Extract multiparameter pooling
      multiparm[[k]] <-
        p.pool
      names(multiparm)[[k]] <-
        paste("Step", k)

      # Clean variable names for selection
      P_temp <-
        clean_P(row.names(p.pool))
      p.pool_temp <-
        p.pool
      row.names(p.pool_temp) <-
        P_temp
      # detect interaction terms and exclude variables
      # that are part of interaction
      if(any(grepl("[*]", P))){
        P_int_id <-
          P_temp[grep("[*]", P_temp)]
        # Detect main effects of interaction terms
        P_main <-
          unique(unlist(str_split(P_int_id, "[*]")))
        # Exclude variables to keep and main effect from selection
        remove_id <-
          !row.names(p.pool_temp) %in% unique(c(P_main))
        p.pool <-
          p.pool[remove_id, , FALSE]
        p.pool_temp <-
          p.pool_temp[remove_id, , FALSE]
      }
      if(!is_empty(keep.P)){
        remove_P_keep <-
          !row.names(p.pool_temp) %in% keep.P
        p.pool <-
          p.pool[remove_P_keep, , FALSE]
      }


      if(p.crit==1){
        break()
      }

      if(nrow(p.pool)==0){
        break()
      }

      # Select variables
      P_temp <-
        row.names(p.pool)
      P_excl <-
        which(p.pool[, 1] == max(p.pool[, 1]))
      if(length(P_excl) > 1) {
        P_excl <-
          P_excl[1]
      }

      if (p.pool[, 1][P_excl] < p.crit) {
        message("\n", "Selection correctly terminated, ",
                "\n", "No more variables removed from the model", "\n")
        (break)()
      }

      P_out <-
        P_temp[P_excl]

      if(p.pool[, 1][P_excl] > p.crit) {
        message("Removed at Step ", k,
                " is - ", P_out)

      }

      P_rm_step[[k]] <-
        P_out
      # Variable excluded on each step
      P <-
        P[!P %in% P_out]

      if(is_empty(P)){
        fm_step[[k+1]] <- as.formula(paste(Y, 1))
        names(fm_step)[[k+1]] <- paste("Step", k+1)
        multiparm[[k+1]] <- 0
        names(multiparm)[[k+1]] <- paste("Step", k+1)
        RR.model[[k+1]] <- print(glm(fm_step[[k+1]], data = data, family = binomial))
        names(RR.model)[[k+1]] <- paste("Step", k+1)
        break()
      }
}
# End k loop
##############################################################

P_rm_step_final <-
  P_rm_step
if(is_empty(P)) {
  P_rm_step <-
    lapply(P_rm_step, clean_P)
} else {
  P_rm_step <-
    lapply(P_rm_step[-k], clean_P)
}

# Extract selected models
outOrder_step <-
  P_orig_temp
if(!is_empty(P_rm_step)){
  P_remove <-
    data.frame(do.call("rbind",
                       lapply(P_rm_step, function(x) {
                         outOrder_step %in% x })))
  names(P_remove) <-
    P_orig
  row.names(P_remove) <-
    paste("Step", 1:length(P_rm_step))
  if(nrow(P_remove)!=1) {
    P_remove <-
      apply(P_remove, 2, function(x) ifelse(x, 1, 0))
    P_remove_final <-
      colSums(P_remove)
    P_remove <-
      rbind(P_remove, P_remove_final)
    row.names(P_remove)[nrow(P_remove)] <-
      "Removed"
  } else {
    P_remove <-
      matrix(apply(P_remove, 2, function(x) ifelse(x, 1, 0)), 1, length(P_orig))
    dimnames(P_remove) <-
      list("Removed", P_orig)
  }
} else {
  P_remove <-
    matrix(rep(0, length(P_orig)), 1, length(P_orig))
  dimnames(P_remove) <-
    list("Removed", P_orig)
}

P_included <-
  as_tibble(names(P_remove[nrow(P_remove), ][P_remove[nrow(P_remove), ] == 0] ))
predictors_final <-
  names(P_remove[nrow(P_remove), ][P_remove[nrow(P_remove), ] == 0])
if(length(P_orig)==1 & !is_empty(P))
  P_included <- predictors_final <- P

if(is_empty(P) | p.crit==1){
  RR_model_step <-
    RR.model
  multiparm_step <-
    multiparm
  fm_step_total <-
    fm_step
  if(p.crit!=1) names(RR_model_step) <- names(multiparm_step) <- names(fm_step_total) <-
    paste("Step", 1:(k+1), "- removal -", c(unlist(P_rm_step_final), "ended"))
  if(p.crit==1)  names(RR_model_step) <- names(multiparm_step) <- names(fm_step_total) <-
    paste("Step", 1, "- no variables removed -")
  RR_model_final <-
    RR_model_step[k+1]
  multiparm_final <-
    multiparm_step[k+1]
  fm_step_final <-
    fm_step_total[k+1]
}
if(!is_empty(P) & p.crit !=1){
  RR_model_step <-
    RR.model
  multiparm_step <-
    multiparm
  fm_step_total <-
    fm_step
  names(RR_model_step) <-
    names(multiparm_step) <- names(fm_step_total) <-
    paste("Step", 1:k, "- removal -", c(unlist(P_rm_step_final), "ended"))
  RR_model_final <-
    RR.model[k]
  multiparm_final <-
    multiparm[k]
  fm_step_final <-
    fm_step[k]
}

Y_initial <-
  c(paste(Outcome, paste("~")))
formula_initial <-
  as.formula(paste(Y_initial, paste(P_orig, collapse = "+")))

pobjbw <-
  list(data = data, RR_model_final = RR_model_final, RR_model = RR_model_step, 
       multiparm_final = multiparm_final,
       multiparm = multiparm_step, formula_step = fm_step_total, 
       formula_final = fm_step_final, formula_initial = formula_initial,
       predictors_in = P_included, predictors_out = P_remove, Outcome = Outcome,
       p.crit = p.crit, call = call, model_type = "binomial",
       predictors_final = predictors_final, predictors_initial = P_orig)
class(pobjbw) <- "smods"
return(pobjbw)
}
