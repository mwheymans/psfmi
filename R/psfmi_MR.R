#' Meng & Rubin pooling method called by psfmi_lr
#'
#' \code{psfmi_MR} Function to pool using Meng & Rubin pooling method
#'
#' @param data Data frame or data matrix with stacked multiple imputed datasets.
#'   The original dataset that contains missing values must be excluded from the dataset.
#' @param nimp A numerical scalar. Number of imputed datasets. Default is 5.
#' @param impvar A character vector. Name of the variable that distinguishes the imputed datasets.
#' @param Outcome Character vector containing the name of the outcome variable.
#' @param P Character vector with the names of the predictor variables.
#'   At least one predictor variable has to be defined.
#' @param p.crit A numerical scalar. P-value selection criterium.
#' @param print.method logical vector. If TRUE full matrix with p-values of all variables according to
#'   chosen method (under method) is shown. If FALSE (default) p-value for categorical variables according
#'   to method are shown and for continuous and dichotomous predictors Rubin’s Rules are used
#'
#' @export
psfmi_MR <-
  function(data, nimp, impvar, Outcome, P, p.crit, print.method)
{

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

    pool.p.val <- matrix(0, length(P), 1)

    for (j in 1:length(P)) {
      cov.nam0 <- P[-j]
      if (length(P) == 1) {
        cov.nam0 <- "1"
      }
      cov.nam.compare <- P[j]
      Y <- c(paste(Outcome, paste("~")))
      form1 <- as.formula(paste(Y, paste(P,
                                            collapse = "+")))
      form0 <- as.formula(paste(Y, paste(cov.nam0, collapse = "+")))
      m1 <- m0 <- nimp
      coef.fit1 <- se.fit1 <- coef.fit0 <- se.fit0 <- list()
      for (i in 1:nimp) {
        dataset <- data[data[impvar] == i, ]
        fit.1 <- glm(form1, data = dataset, family = binomial)
        fit.0 <- glm(form0, data = dataset, family = binomial)
        coef.fit1[[i]] <- summary(fit.1)[[12]][, 1]
        se.fit1[[i]] <- summary(fit.1)[[12]][, 2]
        coef.fit0[[i]] <- summary(fit.0)[[12]][, 1]
        if (length(coef.fit0[[i]]==1)) names(coef.fit0[[i]]) <- "intercept"
        se.fit0[[i]] <- summary(fit.0)[[12]][, 2]
        if (length(se.fit0[[i]]==1)) names(se.fit0[[i]]) <- "intercept"
      }

      coef.fit1.qhat <- do.call("rbind", coef.fit1)
      coef.fit0.qhat <- do.call("rbind", coef.fit0)

      out.res1 <- norm::mi.inference(coef.fit1, se.fit1, 0.95)
      model.res1 <- do.call("cbind", out.res1)
      model.res1 <- model.res1[, -c(3, 7, 8)]
      OR <- exp(model.res1[, 1])
      lower.EXP <- exp(model.res1[, 4])
      upper.EXP <- exp(model.res1[, 5])
      model.res1 <- cbind(model.res1, OR, lower.EXP, upper.EXP)
      model.res1 <- round(model.res1, 4)
      out.res0 <- norm::mi.inference(coef.fit0, se.fit0, 0.95)

      dimQ1 <- length(out.res1$est)
      dimQ2 <- dimQ1 - length(out.res0$est)
      formula1 <- formula(fit.1)
      formula0 <- formula(fit.0)
      devM <- devL <- 0
      for (i in (1:nimp)) {
        dataset <- data[data[impvar] == i, ]
        devL <- devL + LLlogistic(formula1, data = dataset,
          out.res1$est) - LLlogistic(formula0, data = dataset, out.res0$est)
        devM <- devM + LLlogistic(formula1, data = dataset,
          coef.fit1.qhat[i, ]) - LLlogistic(formula0,
          data = dataset, coef.fit0.qhat[i, ])
      }
      devL <- devL/nimp
      devM <- devM/nimp
      rm <- ((nimp + 1)/(dimQ2 * (nimp - 1))) * (devM - devL)
      Dm <- devL/(dimQ2 * (1 + rm))
      v <- dimQ2 * (nimp - 1)
      if (v > 4)
        w <- 4 + (v - 4) * ((1 + (1 - 2/v) * (1/rm))^2)
      else w <- v * (1 + 1/dimQ2) * ((1 + 1/rm)^2)/2

      pool.p.val[j, ] <- round(1 - pf(Dm, dimQ2, w), 5)
      pool.MR <- pool.p.val

      dimnames(pool.MR) <- list(c(P), "Pooled p-value")
      pool.MR <- data.frame(pool.MR)
      names(pool.MR) <- "MR p-values"
    }

    cat("\n", "Pooled model (Rubin's Rules)", "\n")
    model.res <- model.res1
    print(model.res)
    if(print.method){
      cat("\n", "MR Pooled p-values", "\n")
      print(pool.MR)
    }
    id.p.value.RR.f <- grep("factor", row.names(model.res))
    id.p.RR.spl <- grep("rcs", row.names(model.res))
    res.RR <- model.res[-c(1, id.p.value.RR.f, id.p.RR.spl), 3]
    names(res.RR) <- row.names(model.res)[-c(1, id.p.value.RR.f, id.p.RR.spl)]
    pool.MR[names(res.RR), 1] <- res.RR
    cat("\n")
    names(pool.MR) <- "Mixed p-values (MR & RR)"
    print(pool.MR)
}
