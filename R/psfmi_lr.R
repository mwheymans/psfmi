#' Pooling and Predictor selection function for backward or forward selection of
#' Logistic regression models in multiply imputed data.
#'
#' \code{psfmi_lr} Pooling and backward or forward selection of Logistic regression
#' prediction models in multiply imputed data using selection methods RR, D1, D2, D3 and MPR.
#'
#' @param data Data frame with stacked multiple imputed datasets.
#'   The original dataset that contains missing values must be excluded from the
#'   dataset. The imputed datasets must be distinguished by an imputation variable,
#'   specified under impvar, and starting by 1.
#' @param formula A formula object to specify the model as normally used by glm.
#'   See under "Details" and "Examples" how these can be specified.
#' @param nimp A numerical scalar. Number of imputed datasets. Default is 5.
#' @param impvar A character vector. Name of the variable that distinguishes the
#' imputed datasets.
#' @param Outcome Character vector containing the name of the outcome variable.
#' @param predictors Character vector with the names of the predictor variables.
#'   At least one predictor variable has to be defined. Give predictors unique names
#'   and do not use predictor name combinations with numbers as, age2, gnder10, etc.
#' @param cat.predictors A single string or a vector of strings to define the
#' categorical variables. Default is NULL categorical predictors.
#' @param spline.predictors A single string or a vector of strings to define the
#' (restricted cubic) spline variables. Default is NULL spline predictors. See details.
#' @param int.predictors A single string or a vector of strings with the names of the variables that form
#'   an interaction pair, separated by a “:” symbol.
#' @param keep.predictors A single string or a vector of strings including the variables that are forced
#'   in the model during predictor selection. All type of variables are allowed.
#' @param nknots A numerical vector that defines the number of knots for each spline predictor separately.
#' @param p.crit A numerical scalar. P-value selection criterium. A value of 1
#'   provides the pooled model without selection.
#' @param method A character vector to indicate the pooling method for p-values to pool the
#'   total model or used during predictor selection. This can be "RR", D1", "D2", "D3" or "MPR".
#'   See details for more information. Default is "RR".
#' @param direction The direction of predictor selection, "BW" means backward selection and "FW"
#'   means forward selection.

#' @details The basic pooling procedure to derive pooled coefficients, standard errors, 95
#'  confidence intervals and p-values is Rubin's Rules (RR). However, RR is only possible when
#'  the model included continuous or dichotomous variables. Specific procedures are
#'  available when the model also included categorical (> 2 categories) or restricted cubic spline
#'  variables. These pooling methods are: “D1” is pooling of the total covariance matrix,
#'  ”D2” is pooling of Chi-square values, “D3” is pooling Likelihood ratio statistics
#'  (method of Meng and Rubin) and “MPR” is pooling of median p-values (MPR rule).
#'  Spline regression coefficients are defined by using the rcs function for restricted cubic
#'  splines of the rms package. A minimum number of 3 knots as defined under knots is required.
#'
#'  A typical formula object has the form \code{Outcome ~ terms}. Categorical variables has to
#'  be defined as \code{Outcome ~ factor(variable)}, restricted cubic spline variables as
#'  \code{Outcome ~ rcs(variable, 3)}. Interaction terms can be defined as
#'  \code{Outcome ~ variable1*variable2} or \code{Outcome ~ variable1 + variable2 + variable1:variable2}.
#'  All variables in the terms part have to be separated by a "+".
#'
#'@return An object of class \code{pmods} (multiply imputed models) from
#'  which the following objects can be extracted: imputed datasets as \code{data}, final selected
#'  pooled model as \code{RR_model_final}, pooled model at each selection step \code{RR_model},
#'  pooled p-values at final step according to pooling method as \code{multiparm_final}, and
#'  at each step as \code{multiparm}, or \code{multiparm_out} (only when direction = "FW"),
#'  formula object at final step as \code{fm_step_final}, and at each step as \code{fm_step},
#'  predictors included at each selection step as \code{predictors_in}, predictors excluded
#'  at each step as \code{predictors_out}, and name of variable to distinguish imputed datasets
#'  as \code{impvar}, \code{nimp}, \code{Outcome}, \code{method}, \code{p.crit}, \code{call},
#'  \code{model_type}, direction of selection as \code{direction}, \code{predictors_final} for
#'  names of predictors in final selection step and \code{predictors_initial} for names of
#'  predictors in start model.
#'
#' @references Eekhout I, van de Wiel MA, Heymans MW. Methods for significance testing of categorical
#'   covariates in logistic regression models after multiple imputation: power and applicability
#'   analysis. BMC Med Res Methodol. 2017;17(1):129.
#' @references Enders CK (2010). Applied missing data analysis. New York: The Guilford Press.
#' @references Meng X-L, Rubin DB. Performing likelihood ratio tests with multiply-imputed data sets.
#'   Biometrika.1992;79:103-11.
#' @references van de Wiel MA, Berkhof J, van Wieringen WN. Testing the prediction error difference between
#'   2 predictors. Biostatistics. 2009;10:550-60.
#' @references Marshall A, Altman DG, Holder RL, Royston P. Combining estimates of interest in prognostic
#'   modelling studies after multiple imputation: current practice and guidelines. BMC Med Res Methodol.
#'   2009;9:57.
#' @references Van Buuren S. (2018). Flexible Imputation of Missing Data. 2nd Edition. Chapman & Hall/CRC
#'   Interdisciplinary Statistics. Boca Raton.
#' @references EW. Steyerberg (2019). Clinical Prediction MOdels. A Practical Approach
#'  to Development, Validation, and Updating (2nd edition). Springer Nature Switzerland AG.
#'
#' @references http://missingdatasolutions.rbind.io/
#' 
#' @author Martijn Heymans, 2020
#'
#' @examples
#'   pool_lr <- psfmi_lr( data=lbpmilr, formula = Chronic ~ Pain + 
#'   factor(Satisfaction) + rcs(Tampascale,3) + Radiation + 
#'   Radiation*factor(Satisfaction) + Age + Duration + BMI,
#'   p.crit = 0.05, direction="FW", nimp=5, impvar="Impnr", 
#'   keep.predictors = c("Radiation*factor(Satisfaction)", "Age"), method="D1")
#'   
#'   pool_lr$RR_model_final
#' 
#' @export
psfmi_lr <- function(data,
                     formula = NULL,
                     nimp=5,
                     impvar=NULL,
                     Outcome = NULL,
                     predictors=NULL,
                     cat.predictors=NULL,
                     spline.predictors=NULL,
                     int.predictors=NULL,
                     keep.predictors=NULL,
                     nknots=NULL,
                     p.crit=1,
                     method="RR",
                     direction=NULL)
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
  if(p.crit!=1){
    if(is_empty(direction))
      stop("Specify FW or BW for forward or backward predictor selection")
  }
  if (!(is.data.frame(data)))
    stop("Data should be a data frame")
  data <- data.frame(as_tibble(data))
  data <- mutate_if(data, is.factor, ~ as.numeric(as.character(.x)))
  if(!all(data[Outcome]==1 | data[Outcome]==0))
    stop("Outcome should be a 0 - 1 variable")
  if ((nvar <- ncol(data)) < 2)
    stop("Data should contain at least two columns")
  if(is_empty(impvar))
    stop("Imputation variable is not defined")
  if(is_empty(method)) method="RR"
  if(all(!is_empty(cat.P) | !is_empty(s.P)) & method=="RR")
    stop("Categorical or spline variables in model, define selection method: D1, D2, D3 or MPR")
  if (order(unique(data[, impvar]))[1] == 0)
    stop("Original dataset should not be included")
  if(is_empty(nimp))
    stop("Number of imputed datasets is not defined, use nimp!")
  if (nimp < 2) {
    stop("\n", "Number of imputed datasets must be > 1", "\n\n")
  }
  if (p.crit > 1)
    stop("\n", "P-value criterium > 1", "\n")
  if (any(nknots<3))
    stop("\n", "Number of knots must be > 2", "\n")
  if (length(nknots) != length(s.P))
    stop("\n", "Number of knots not specified for every spline variable", "\n")
  if (!is_empty(cat.P)) {
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
  if (is_empty(P))
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
      if(!is_empty(keep.P)){
        keep.P <-
          gsub(cat.P,
               replacement=paste0("factor(", cat.P, ")"), keep.P)
      }
    } else {
      for(i in 1:length(cat.P)) {
        P <-
          gsub(cat.P[i],
               replacement=paste0("factor(", cat.P[i], ")"), P)
        if(!is_empty(keep.P)){
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
      if(!is_empty(keep.P)){
        keep.P <-
          gsub(s.P,
               replacement=paste0("rcs(", s.P, ",", nknots, ")"), keep.P)
      }
    } else {
      for(i in 1:length(s.P)) {
        P <- gsub(s.P[i],
                  replacement=paste0("rcs(", s.P[i], ",", nknots[i], ")"), P)
        if(!is_empty(keep.P)){
          keep.P <-
            gsub(s.P[i],
                 replacement=paste0("rcs(", s.P[i], ",", nknots[i], ")"), keep.P)
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
  
  if(any(!keep.P %in% P))
    stop("\n", "Variables to keep not defined as Predictor", "\n\n")
  
  if(p.crit==1){
    pobjpool <-
      psfmi_lr_bw(data = data, nimp=nimp, impvar = impvar, Outcome= Outcome,
                  P = P, p.crit = p.crit, method = method, keep.P = keep.P)
    class(pobjpool) <-
      "pmods"
    return(pobjpool)
  }
  if(direction=="FW"){
    pobjfw <-
      psfmi_lr_fw(data = data, nimp = nimp, impvar = impvar, Outcome = Outcome, p.crit = p.crit,
                  P = P, keep.P = keep.P, method = method)
    class(pobjfw) <-
      "pmods"
    return(pobjfw)
  }
  if(direction=="BW"){
    pobjbw <-
      psfmi_lr_bw(data = data, nimp=nimp, impvar = impvar, Outcome= Outcome,
                  P = P, p.crit = p.crit, method = method, keep.P = keep.P)
    class(pobjbw) <-
      "pmods"
    return(pobjbw)
  }
}