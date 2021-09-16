#' Backward selection of Logistic regression models in multiply imputed data.
#'
#' \code{psfmi_lr_bw} Backward selection of Logistic regression
#' models in multiply imputed data using selection methods RR, D1, D2, D3 and MPR.
#' Function is called by \code{psfmi_lr}. 
#'
#' @param data Data frame with stacked multiple imputed datasets.
#'   The original dataset that contains missing values must be excluded from the
#'   dataset. The imputed datasets must be distinguished by an imputation variable,
#'   specified under impvar, and starting by 1.
#' @param nimp A numerical scalar. Number of imputed datasets. Default is 5.
#' @param impvar A character vector. Name of the variable that distinguishes the
#' imputed datasets.
#' @param Outcome Character vector containing the name of the outcome variable.
#' @param P Character vector with the names of the predictor variables.
#'   At least one predictor variable has to be defined. Give predictors unique names
#'   and do not use predictor name combinations with numbers as, age2, BMI10, etc.
#' @param p.crit A numerical scalar. P-value selection criterium. A value of 1
#'   provides the pooled model without selection.
#' @param method A character vector to indicate the pooling method for p-values to pool the
#'   total model or used during predictor selection. This can be "RR", D1", "D2", "D3" or "MPR".
#'   See details for more information. Default is "RR".
#' @param keep.P A single string or a vector of strings including the variables that are forced
#'   in the model during predictor selection. All type of variables are allowed.
#'   
#' @author Martijn Heymans, 2020
#' @keywords internal
#'  
#' @export
psfmi_lr_bw <- function(data, nimp, impvar, Outcome, P, p.crit, method, keep.P)
{
  call <- match.call()

  RR.model <- P_rm_step <- fm_step <- imp.dt <- multiparm <- list()

  P_orig <-
    P
  keep.P <-
    sapply(as.list(keep.P), clean_P)
  P_orig_temp <-
    clean_P(P)

  if(!is_empty(keep.P))
    if(length(P_orig)==1)
      if(P_orig_temp == keep.P)
        stop("\n", "No need to define keep.predictors. Exclude keep.predictors and set p.crit = 1","\n")

  # Loop k, to pool models in multiply imputed datasets
  for (k in 1:(length(P)+1)) {

    # set regression formula fm
    Y <-
      c(paste(Outcome, paste("~")))
    fm <-
      as.formula(paste(Y, paste(P, collapse = "+")))

    # Extract df of freedom for MPR
    if(method=="MPR" | method=="RR"){

      # Extract LR and P value of model for MPR
      chi.LR <-
        data.frame(matrix(0, length(attr(terms(fm), "term.labels")), nimp))
      chi.p <-
        data.frame(matrix(0, length(attr(terms(fm), "term.labels")), nimp))

      fit <- list()
      for (i in 1:nimp) {
        imp.dt[[i]] <- data[data[impvar] == i, ]
        fit[[i]] <- glm(fm, data = imp.dt[[i]], family = binomial)
        if(length(attr(terms(fm), "term.labels")) == 1){
          chi.LR[, i] <- car::Anova(fit[[i]])[1]
          chi.p[, i] <- car::Anova(fit[[i]])[3]
        } else {
          chi.LR[, i] <- car::Anova(fit[[i]])[, 1]
          chi.p[, i] <- car::Anova(fit[[i]])[, 3]
        }
      }

      # Rubin's Rules
      out.res <-
        summary(pool(fit))
      OR <-
        exp(out.res$estimate)
      lower.EXP <-
        exp(out.res$estimate - (qt(0.975, out.res$df)*out.res$std.error))
      upper.EXP <-
        exp(out.res$estimate + (qt(0.975, out.res$df)*out.res$std.error))
      model.res <-
        data.frame(cbind(out.res, OR, lower.EXP, upper.EXP))
      RR.model[[k]] <-
        model.res
      names(RR.model)[[k]] <-
        paste("Step", k)#
    }

    # D1 and D2 pooling methods
    if(method=="D1" | method == "D2" | method=="D3" | method=="D4") {

      pool.p.val <-
        matrix(0, length(P), 2)
      P_test <-
        clean_P(P)

      for (j in 1:length(P)) {
        cov.nam0 <-
          P[-j]
        if (length(P) == 1) {
          cov.nam0 <-
            "1"
        }
        Y <-
          c(paste(Outcome, paste("~")))
        form1 <-
          as.formula(paste(Y, paste(P, collapse = "+")))
        form0 <-
          as.formula(paste(Y, paste(cov.nam0, collapse = "+")))

        if(any(grepl(P_test[j], P_test[-j]))){
          cov.nam0 <-
            P[-grep(P_test[j], P_test)]
          form0 <-
            as.formula(paste(Y, paste(c(cov.nam0), collapse = "+")))
        }
        if(method=="D4"){
          data <-
            filter(data, data[impvar] <= nimp)
          imp_list <-
            data %>% group_split(data[, impvar], .keep = FALSE) %>%
            mitools::imputationList(imp_list)
          
          fit0 <-
            with(data=imp_list, expr= glm(as.formula(paste(Y,
                       paste(cov.nam0, collapse = "+"))), family = binomial))
          fit1 <-
            with(data=imp_list, expr= glm(as.formula(paste(Y,
                       paste(P, collapse = "+"))), family = binomial))
          
          out.res1 <-
            summary(pool(fit1))
          OR <-
            exp(out.res1$estimate)
          lower.EXP <-
            exp(out.res1$estimate - (qt(0.975, out.res1$df)*out.res1$std.error))
          upper.EXP <-
            exp(out.res1$estimate + (qt(0.975, out.res1$df)*out.res1$std.error))
          model.res1 <-
            data.frame(cbind(out.res1, OR, lower.EXP, upper.EXP))
          RR.model[[k]] <-
            model.res1
          names(RR.model)[[k]] <-
            paste("Step", k)
          
          res_D4 <- 
            pool_D4(data=data, fm0=form0, fm1=form1, nimp=nimp,
                            impvar=impvar, robust=TRUE, model_type="binomial")
          pvalue <-
            res_D4$pval
          fstat <-
            res_D4$F
          pool.p.val[j, ] <-
            c(pvalue, fstat)
        }
        if(method=="D3"){
          data <-
            filter(data, data[impvar] <= nimp)
          imp_list <-
            data %>% group_split(data[, impvar], .keep = FALSE) %>%
            mitools::imputationList(imp_list)

          fit0 <-
            with(data=imp_list, expr= glm(as.formula(paste(Y,
                   paste(cov.nam0, collapse = "+"))), family = binomial))
          fit1 <-
            with(data=imp_list, expr= glm(as.formula(paste(Y,
                   paste(P, collapse = "+"))), family = binomial))

          out.res1 <-
            summary(pool(fit1))
          OR <-
            exp(out.res1$estimate)
          lower.EXP <-
            exp(out.res1$estimate - (qt(0.975, out.res1$df)*out.res1$std.error))
          upper.EXP <-
            exp(out.res1$estimate + (qt(0.975, out.res1$df)*out.res1$std.error))
          model.res1 <-
            data.frame(cbind(out.res1, OR, lower.EXP, upper.EXP))
          RR.model[[k]] <-
            model.res1
          names(RR.model)[[k]] <-
            paste("Step", k)

          test_P <-
            mice::D3(fit1, fit0)
          pvalue <-
            test_P$result[4]
          fstat <-
            test_P$result[1]
          pool.p.val[j, ] <-
            c(pvalue, fstat)
        }
        if(method =="D1" | method =="D2"){
          fit1 <- fit0 <- imp.dt <- list()
          for (i in 1:nimp) {
            imp.dt[[i]] <-
              data[data[impvar] == i, ]
            fit1[[i]] <-
              glm(form1, data = imp.dt[[i]], family = binomial)
            fit0[[i]] <-
              glm(form0, data = imp.dt[[i]], family = binomial)
          }
          out.res1 <-
            summary(pool(fit1))
          OR <-
            exp(out.res1$estimate)
          lower.EXP <-
            exp(out.res1$estimate - (qt(0.975, out.res1$df)*out.res1$std.error))
          upper.EXP <-
            exp(out.res1$estimate + (qt(0.975, out.res1$df)*out.res1$std.error))
          model.res1 <-
            data.frame(cbind(out.res1, OR, lower.EXP, upper.EXP))
          RR.model[[k]] <-
            model.res1
          names(RR.model)[[k]] <-
            paste("Step", k)

          test_P <-
            mitml::testModels(fit1, fit0, method = method)
          pvalue <-
            test_P$test[4]
          fstat <-
            test_P$test[1]
          pool.p.val[j, ] <-
            c(pvalue, fstat)
        }
      }
      p.pool <-
        data.frame(pool.p.val)
      row.names(p.pool) <-
        P
      names(p.pool) <-
        c(paste("p-values", method), "F-statistic")
    }

    # MPR Pooling
    if(method=="MPR") {
      p.pool <-
        data.frame(apply(chi.p, 1 , median))
      rownames(p.pool) <-
        P
      names(p.pool) <-
        c("p-value MPR")
    }

    # RR Pooling
    if(method=="RR") {
      p.pool <-
        data.frame(RR.model[[k]][-1, 6],
                   row.names=P)
      names(p.pool) <-
        c("p-value RR")
    }

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

    if(nrow(p.pool)==0)
      break()

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
      fm_step[[k+1]] <-
        as.formula(paste(Y, 1))
      names(fm_step)[[k+1]] <-
        paste("Step", k+1)
      multiparm[[k+1]] <-
        0
      names(multiparm)[[k+1]] <-
        paste("Step", k+1)
      fit <- list()
      for (i in 1:nimp) {
        imp.dt[[i]] <-
          data[data[impvar] == i, ]
        fit[[i]] <-
          glm(fm_step[[k+1]], data = imp.dt[[i]], family = binomial)
      }

      # Rubin's Rules
      out.res <-
        summary(pool(fit))
      OR <-
        exp(out.res$estimate)
      lower.EXP <-
        exp(out.res$estimate - (qt(0.975, out.res$df)*out.res$std.error))
      upper.EXP <-
        exp(out.res$estimate + (qt(0.975, out.res$df)*out.res$std.error))
      model.res <-
        data.frame(cbind(out.res, OR, lower.EXP, upper.EXP))
      RR.model[[k+1]] <-
        model.res
      names(RR.model)[[k+1]] <-
        paste("Step", k+1)
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
    if(p.crit==1) {
      Y_initial <-
        c(paste(Outcome, paste("~")))
      formula_initial <-
        as.formula(paste(Y_initial, paste(P_orig, collapse = "+")))
      fm_step_final <- 
        formula_initial
      RR_model_final <-
        RR_model_step[k]
      multiparm_final <-
        multiparm_step[k]
    }
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

  bw <-
    list(data = data, RR_model = RR_model_step, RR_model_final = RR_model_final, 
         multiparm = multiparm_step, multiparm_final = multiparm_final, 
         formula_step = fm_step_total, formula_final = fm_step_final,
         formula_initial = formula_initial, 
         predictors_in = P_included, predictors_out = P_remove,
         impvar = impvar, nimp = nimp, Outcome = Outcome,
         method = method, p.crit = p.crit, call = call, 
         model_type = "binomial", direction = "BW", 
         predictors_final = predictors_final,
         predictors_initial = P_orig, keep.predictors = keep.P)
  return(bw)
}