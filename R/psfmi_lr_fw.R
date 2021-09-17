#' Forward selection of Logistic regression models in multiply imputed data.
#'
#' \code{psfmi_lr_fw} Forward selection of Logistic regression
#' models across multiply imputed data using selection methods RR, D1, D2, D3 and MPR.
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
#'   in the model during predictor selection. Categorical and interaction variables are allowed.
#'   
#' @author Martijn Heymans, 2020
#' @keywords internal
#'    
#' @export
psfmi_lr_fw <- function(data, nimp, impvar, Outcome, P, p.crit, method, keep.P)
{

  call <- match.call()

  P_each_step <- fm_step <- fm_total <- RR_model_total <-
    RR_model_select <- imp.dt <- multiparm_step <- multiparm_end <- list()
  
  P_select <- 0
  P_orig <- P
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
    keep.P <- keep.P_temp
  }
  
  # Start J loop, to build up models,
  # variable by variable
  for(j in 1:length(P))
  {
    
    P_D1 <- P_D2 <- P_D3 <- P_D4 <- P_MPR <- P_RR <-
      RR.model <- fm_step <- as.list(rep(0, length(P)))
    
    # Loop k, to pool models in multiply imputed datasets
    for (k in 1:length(P)) {
      
      # set regression formula fm
      Y <- c(paste(Outcome, paste("~")))
      fm <- terms(as.formula(paste0(Y, paste0(c(P[k], keep.P), collapse = "+"))))
      
      if(P_select!=0){
        fm <- terms(update.formula(fm, paste0("~. +",
                                              paste0(paste0(P_each_step, collapse = "+")))))
      }
      
      # Extract df of freedom for MPR
      if(method=="MPR" | method=="RR"){
        
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
        out.res <- summary(pool(fit))
        OR <- exp(out.res$estimate)
        lower.EXP <- exp(out.res$estimate - (qt(0.975, out.res$df)*out.res$std.error))
        upper.EXP <- exp(out.res$estimate + (qt(0.975, out.res$df)*out.res$std.error))
        model.res <- data.frame(cbind(out.res, OR, lower.EXP, upper.EXP))
        RR.model[[k]] <- model.res
      }
      
      # D1 and D2 pooling methods
      if(method=="D1" | method == "D2" | method == "D4") {
        
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
            pool_D4(data=data, fm0=f0, fm1=fm, nimp=nimp,
                    impvar=impvar, robust=TRUE, model_type="binomial")
          pvalue <-
            res_D4$pval
          fstat <-
            res_D4$F
          pool.multiparm <- data.frame(matrix(c(pvalue, fstat), length(P[k]), 2))
          row.names(pool.multiparm) <- P[k]
          names(pool.multiparm) <- c("p-values", "F-statistic")
          pool.multiparm
        }
        if(method=="D1" | method == "D2"){
          fit1 <- fit0 <- list()
          for (i in 1:nimp) {
            imp.dt[[i]] <- data[data[impvar] == i, ]
            fit1[[i]] <- glm(fm, data = imp.dt[[i]], family = binomial)
            fit0[[i]] <- glm(f0, data = imp.dt[[i]], family = binomial)
          }
          
          out.res <- summary(pool(fit1))
          OR <- exp(out.res$estimate)
          lower.EXP <- exp(out.res$estimate - (qt(0.975, out.res$df)*out.res$std.error))
          upper.EXP <- exp(out.res$estimate + (qt(0.975, out.res$df)*out.res$std.error))
          model.res <- data.frame(cbind(out.res, OR, lower.EXP, upper.EXP))
          RR.model[[k]] <- model.res
          names(RR.model)[k] <- paste("Step", j)
          if(P_select==0) names(RR.model)[k] <- paste("Step", 1)
          
          tmr <- mitml::testModels(fit1, fit0, method = method)
          
          pvalue <- tmr$test[4]
          fstat <- tmr$test[1]
          pool.multiparm <- data.frame(matrix(c(pvalue, fstat), length(P[k]), 2))
          row.names(pool.multiparm) <- P[k]
          names(pool.multiparm) <- c("p-values", "F-statistic")
          pool.multiparm
        }
        
        if(method=="D1") P_D1[[k]] <- pool.multiparm
        if(method=="D2") P_D2[[k]] <- pool.multiparm
        if(method=="D4") P_D4[[k]] <- pool.multiparm
        
        # Set f0 to original, before testing interactions
        if(P_select==0){
          cov.nam0 <- "1"
          if(!is_empty(keep.P)) cov.nam0 <- keep.P
          f0 <- as.formula(paste(Y, paste(c(cov.nam0), collapse = "+")))
        }
        if(P_select!=0)
          f0 <- f0_orig
        
      }
      # D3 pooling
      if(method=="D3") {
        
        LLlogistic <-
          function(formula, data, coefs) {
            logistic <- function(mu) exp(mu)/(1 + exp(mu))
            Xb <- model.matrix(formula, data) %*% coefs
            y <- model.frame(formula, data)[1][, 1]
            p <- logistic(Xb)
            y <- (y - min(y))/(max(y) - min(y))
            term1 <- term2 <- rep(0, length(y))
            term1[y != 0] <- y[y != 0] * log(y[y != 0]/p[y != 0])
            term2[y == 0] <- (1 - y[y == 0]) * log((1 - y[y == 0])/(1 - p[y == 0]))
            return(-(2 * sum(term1 + term2)))
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
        
        m1 <- m0 <- nimp
        fit1 <- fit0 <- coef.fit1 <- se.fit1 <- coef.fit0 <- se.fit0 <- list()
        for (i in 1:nimp) {
          dataset <- data[data[impvar] == i, ]
          fit1[[i]] <- glm(fm, data = dataset, family = binomial)
          fit0[[i]] <- glm(f0, data = dataset, family = binomial)
          coef.fit1[[i]] <- summary(fit1[[i]])[[12]][, 1]
          se.fit1[[i]] <- summary(fit1[[i]])[[12]][, 2]
          coef.fit0[[i]] <- summary(fit0[[i]])[[12]][, 1]
          if (length(coef.fit0[[i]])==1) names(coef.fit0[[i]]) <- "intercept"
          se.fit0[[i]] <- summary(fit0[[i]])[[12]][, 2]
          if (length(se.fit0[[i]])==1) names(se.fit0[[i]]) <- "intercept"
        }
        
        coef.fit1.qhat <- do.call("rbind", coef.fit1)
        coef.fit0.qhat <- do.call("rbind", coef.fit0)
        
        out.res1 <- summary(pool(fit1))
        OR <- exp(out.res1$estimate)
        lower.EXP <- exp(out.res1$estimate - (qt(0.975, out.res1$df)*out.res1$std.error))
        upper.EXP <- exp(out.res1$estimate + (qt(0.975, out.res1$df)*out.res1$std.error))
        model.res1 <- data.frame(cbind(out.res1, OR, lower.EXP, upper.EXP))
        RR.model[[k]] <- model.res1
        names(RR.model)[k] <- paste("Step", j)
        if(P_select==0) names(RR.model)[k] <- paste("Step", 1)
        
        out.res0 <- summary(pool(fit0))
        dimQ1 <- length(out.res1$est)
        dimQ2 <- dimQ1 - length(out.res0$est)
        formula1 <- formula(fit1[[1]])
        formula0 <- formula(fit0[[1]])
        devM <- devL <- 0
        for (i in (1:nimp)) {
          dataset <- data[data[impvar] == i, ]
          devL <- devL + LLlogistic(formula1, data = dataset,
                                    out.res1$est) - LLlogistic(formula0,
                                                               data = dataset, out.res0$est)
          devM <- devM + LLlogistic(formula1, data = dataset,
                                    coef.fit1.qhat[i, ]) - LLlogistic(formula0,
                                                                      data = dataset, coef.fit0.qhat[i, ])
        }
        devL <- devL/nimp
        devM <- devM/nimp
        rm <- ((nimp + 1)/(dimQ2 * (nimp - 1))) * (devM - devL)
        Dm <- devL/(dimQ2 * (1 + rm))
        v <- dimQ2 * (nimp - 1)
        if (v > 4) {
          w <- 4 + (v - 4) * ((1 + (1 - 2/v) * (1/rm))^2)
        } else {
          w <- v * (1 + 1/dimQ2) * ((1 + 1/rm)^2)/2
        }
        pool.p.val <- round(1 - pf(Dm, dimQ2, w), 5)
        pool.multiparm <- data.frame(matrix(c(pool.p.val, Dm), length(P[k]), 2))
        row.names(pool.multiparm) <- P[k]
        names(pool.multiparm) <- c("p-values", "D3 statistic")
        pool.multiparm
        
        P_D3[[k]] <- pool.multiparm
        
        # Set f0 to original, before testing interactions
        if(P_select==0){
          cov.nam0 <- "1"
          if(!is_empty(keep.P)) cov.nam0 <- keep.P
          f0 <- as.formula(paste(Y, paste(c(cov.nam0), collapse = "+")))
        }
        if(P_select!=0)
          f0 <- f0_orig
      }
      
      # MPR Pooling
      if(method=="MPR") {
        med.pvalue <- data.frame(apply(chi.p, 1 , median))
        rownames(med.pvalue) <- clean_P(attr(fm, "term.labels")) %>%
          str_replace(":", "*")
        names(med.pvalue) <- c("p-value MPR")
        P_MPR[[k]] <- med.pvalue
      }
      
      if(method=="RR") {
        RR <- data.frame(RR.model[[k]])[, c(1,6)]
        names(RR)[2] <- c("p-value RR")
        P_RR[[k]] <- RR
      }
      
      # Extract regression formula's
      fm_step[[k]] <- paste(Y, paste(attr(fm, "term.labels"), collapse = " + "))
      names(fm_step)[k] <- paste("Test - ", P[k])
    }
    # End k loop
    ##############################################################
    
    fm_total[[j]] <- fm_step
    RR_model_total[[j]] <- RR.model
    
    # p.pool for RR
    if(method=="RR"){
      P_RR_id <- P
      P_RR_id <- P_RR_id %>% str_replace("[*]", ":")
      p.pool <- data.frame("Pvalue"=do.call("rbind", purrr::pmap(list(x=P_RR, y=as.list(P_RR_id)),
                                                                 function(x, y) { x[x[, "term"] == y, -1] }) ))
      row.names(p.pool) <- P_RR_id
      names(p.pool) <- paste("p-value", method)
    }
    
    # p.pool for multiparameer pooling D1, D2, D3
    if(method=="D1" | method == "D2" | method == "D3" | method == "D4"){
      p.pool <- data.frame(do.call("rbind", get(paste("P", method, sep="_"))))
      rownames_temp <- row.names(p.pool)
      p.pool <- data.frame(p.pool[, 1])
      row.names(p.pool) <- clean_P(rownames_temp)
      names(p.pool) <- paste("p-value", method)
    }
    
    # p.pool for MPR
    if(method=="MPR"){
      P_MPR_id <- P
      P_MPR_id <- clean_P(P_MPR_id)
      p.pool <- do.call("rbind", purrr::pmap(list(x=P_MPR, y=as.list(P_MPR_id)),
                                             function(x, y) {
                                               x <- data.frame(x[row.names(x)==y,])
                                               row.names(x) <- y
                                               names(x) <- "P-value"
                                               return(x)
                                             }) )
    }
    
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
      P_each_step <- c(P_each_step[-j])#, keep.P)
      if(is_empty(P_each_step)){
        fm_step <- as.formula(paste(Y, 1))
        if(!is_empty(keep.P)){
          fm_step <- as.formula(paste(Y, paste(keep.P, collapse = "+")))
          P_each_step <- c(keep.P)
        }
        fit <- list()
        for (i in 1:nimp) {
          imp.dt[[i]] <- data[data[impvar] == i, ]
          fit[[i]] <- glm(fm_step, data = imp.dt[[i]], family = binomial)
        }
        # Rubin's Rules
        out.res <- summary(pool(fit))
        OR <- exp(out.res$estimate)
        lower.EXP <- exp(out.res$estimate - (qt(0.975, out.res$df)*out.res$std.error))
        upper.EXP <- exp(out.res$estimate + (qt(0.975, out.res$df)*out.res$std.error))
        model.res <- data.frame(cbind(out.res, OR, lower.EXP, upper.EXP))
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
      multiparm_final <- multiparm[j]
      fm_step_final <- fm_total[j]
    }
    if(!is_empty(P)){
      RR_model_final <- RR_model[j-1]
      multiparm_final <- multiparm[j-1]
      fm_step_final <- fm_total[j-1]
      if(j==1 & !is_empty(keep.P)) {
        RR_model_final <- RR_model
        fm_step_final <- fm_total
      }
      names(RR_model_final) <- "Final model"
      multiparm_out <- multiparm_end[j]
      names(multiparm_out) <- "Predictors removed"
    }
  }
  
  
  Y_initial <-
    c(paste(Outcome, paste("~")))
  formula_initial <-
    as.formula(paste(Y_initial, paste(P_orig, collapse = "+")))
  
  fw <- list(data = data, RR_model = RR_model, RR_model_final = RR_model_final,
             multiparm = multiparm, multiparm_final = multiparm_final,
             multiparm_step = multiparm_step, multiparm_out = multiparm_out,
             formula_step = fm_total, formula_final = fm_step_final,
             formula_initial = formula_initial,
             predictors_in = P_select, predictors_out = P_excluded,
             impvar = impvar, nimp = nimp, Outcome = Outcome,
             method = method, p.crit = p.crit, call = call,
             model_type = "binomial", direction = "FW",
             predictors_final = predictors_final, predictors_initial = P_orig,
             keep.predictors = keep.P)
  
  return(fw)
}