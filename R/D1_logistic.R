#' D1 method for Predictor selection called by psfmi_lr
#'
#' \code{D1_logistic} D1 pooling method
#'
#' @param data Data frame or data matrix with stacked multiple imputed datasets.
#'   The original dataset that contains missing values must be excluded from the dataset.
#' @param nimp A numerical scalar. Number of imputed datasets. Default is 5.
#' @param impvar A character vector. Name of the variable that distinguishes the imputed datasets.
#' @param fm regression formula from glm object
#' @param names.var list of predictors included in pooled regression model
#'
#' @export
D1_logistic <-
  function(data, impvar, nimp, fm, names.var)
{

  # Regression equation

  fa <- glm(fm, family=binomial, data=data[data[impvar]==1, ])

  k <- length(coef(fa))
  names <- names(coef(fa))

  qhat <- matrix(NA, nrow = nimp, ncol = k, dimnames = list(1:nimp, names))
  u <- array(NA, dim = c(nimp, k, k), dimnames = list(1:nimp, names, names))

  for (i in 1:nimp)
  {
    # Regression equation
    dataset <- data[data[impvar]==i, ]
    fit <- glm(fm, family=binomial, data=dataset)
    qhat[i, ] <- coef(fit)
    ui <- vcov(fit)
    u[i, , ] <- ui
  }

  # Average coefficient estimate
  qbar <- apply(qhat, 2, mean)

  # Within imputation covariance
  ubar <- apply(u, c(2, 3), mean)
  e <- qhat - matrix(qbar, nrow = nimp, ncol = k, byrow = TRUE)

  # Between imputation covariance
  b <- (t(e) %*% e)/(nimp - 1)

  res <- list()

  for(i in 1:length(names.var))
  {
    nr.var <- length(names.var[[i]])
    r <- ((1 + 1/nimp) * sum(diag(b/ubar)[names.var[[i]]]))/nr.var
    # More stable estimate of total variance
    T.stable <- (1+r)*ubar[names.var[[i]], names.var[[i]]]
    # D1 Statistic
    theta0 <- matrix(0, nrow = nr.var, ncol = 1)
    thetahat <- matrix(qbar[names.var[[i]]] , nrow=nr.var, ncol=1)
    z <- thetahat - theta0
    z

    d1 <- 1/nr.var * (t(z) %*% (solve(T.stable)%*%z))
    d1 <- round(d1, 3)
    t <- nr.var*(nimp-1)
    t
    if(t < 5) {
      v1 <- 1/2 * (nr.var + 1) * (nimp-1) * (1 + 1/r)^2
      v1
      p.value <- round(1-stats::pf(d1, nr.var, v1), 4)
    }
    if (t>4) {
      v2 <- 4 + (t-4) * ((1 + (1/r * (1-2/t)))^2)
      round(v2, 0)
      p.value <- round(1-stats::pf(d1, nr.var, v2), 4)
    }
    res[[i]] <- c("F-statistic"=d1, "P-value"=p.value)
  }
  #print(res)
  res <- matrix(unlist(res), length(names.var), 2, byrow=T)
  dimnames(res) <- list(names.var, c("Chi-sq", "Pooled P-value"))
  res

}
