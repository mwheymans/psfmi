#' Function to pool NRI measures over Multiply Imputed datasets  
#' 
#' \code{pool_reclassification} Function to pool categorical and continuous NRI
#'   and IDI over Multiply Imputed datasets 
#'
#'@param datasets a list of data frames corresponding to the multiply imputed
#'  datasets, within each dataset in the first column the predicted probabilities
#'  of model 1, in the second column those of model 2 and in the third column
#'  the observed outcomes coded as '0'and '1'.
#'@param cutoff cutoff value for the categorical NRI, must lie between 0 and 1.
#'  
#'@details This function is called by the function \code{pool_compare_model}  
#'
#'@author Martijn Heymans, 2020
#'
#'@export
pool_reclassification <-
  function(datasets, cutoff = cutoff) {
    
    if(!inherits(datasets, "list"))
      stop("\n", "Object should be a list of data.frames or matrices
               with predicted probabilities and outcomes corresponding
               to the multiply imputed datasets", "\n")
    
    cutoff <-
      c(0, cutoff, 1)
    reclass <-
      lapply(datasets, function(x) {
        
        pred1 <-
          x[, 1]
        pred2 <-
          x[, 2]
        y <-
          x[, 3]
        
        # NRI categorical
        c1 <-
          cut(pred1, breaks = cutoff, include.lowest = TRUE, right = FALSE)
        c2 <-
          cut(pred2, breaks = cutoff, include.lowest = TRUE, right = FALSE)
        c11 <-
          factor(c1, levels = levels(c1), labels = c(1:length(levels(c1))))
        c22 <-
          factor(c2, levels = levels(c2), labels = c(1:length(levels(c2))))
        
        pred1_cat <-
          as.numeric(c11) * (1/(length(levels(c11))))
        pred2_cat <-
          as.numeric(c22) * (1/(length(levels(c22))))
        
        n <-
          length(y)
        y <-
          as.numeric(y)
        u <-
          sort(unique(y))
        a <-
          y == 1
        b <-
          y == 0
        na <-
          sum(a)
        nb <-
          sum(b)
        
        # nri categorical
        d_cat <- pred2_cat - pred1_cat
        nup.ev_cat <- sum(d_cat[a] > 0)
        pup.ev_cat <- nup.ev_cat/na
        nup.ne_cat <- sum(d_cat[b] > 0)
        pup.ne_cat <- nup.ne_cat/nb
        ndown.ev_cat <- sum(d_cat[a] < 0)
        pdown.ev_cat <- ndown.ev_cat/na
        ndown.ne_cat <- sum(d_cat[b] < 0)
        pdown.ne_cat <- ndown.ne_cat/nb
        nri.ev_cat <- pup.ev_cat - pdown.ev_cat
        v.nri.ev_cat <- (nup.ev_cat + ndown.ev_cat)/(na^2) - ((nup.ev_cat - ndown.ev_cat)^2)/(na^3)
        se.nri.ev_cat <- sqrt(v.nri.ev_cat)
        #z.nri.ev <- nri.ev/se.nri.ev
        nri.ne_cat <- pdown.ne_cat - pup.ne_cat
        v.nri.ne_cat <- (ndown.ne_cat + nup.ne_cat)/(nb^2) - ((ndown.ne_cat - nup.ne_cat)^2)/(nb^3)
        se.nri.ne_cat <- sqrt(v.nri.ne_cat)
        #z.nri.ne <- nri.ne/se.nri.ne
        nri_cat <- pup.ev_cat - pdown.ev_cat - (pup.ne_cat - pdown.ne_cat)
        se.nri_cat <- sqrt(v.nri.ev_cat + v.nri.ne_cat)
        
        # nri continuous and idi
        d <- pred2 - pred1
        nup.ev <- sum(d[a] > 0)
        pup.ev <- nup.ev/na
        nup.ne <- sum(d[b] > 0)
        pup.ne <- nup.ne/nb
        ndown.ev <- sum(d[a] < 0)
        pdown.ev <- ndown.ev/na
        ndown.ne <- sum(d[b] < 0)
        pdown.ne <- ndown.ne/nb
        nri.ev <- pup.ev - pdown.ev
        v.nri.ev <- (nup.ev + ndown.ev)/(na^2) - ((nup.ev - ndown.ev)^2)/(na^3)
        se.nri.ev <- sqrt(v.nri.ev)
        z.nri.ev <- nri.ev/se.nri.ev
        nri.ne <- pdown.ne - pup.ne
        v.nri.ne <- (ndown.ne + nup.ne)/(nb^2) - ((ndown.ne - nup.ne)^2)/(nb^3)
        se.nri.ne <- sqrt(v.nri.ne)
        #z.nri.ne <- nri.ne/se.nri.ne
        nri <- pup.ev - pdown.ev - (pup.ne - pdown.ne)
        se.nri <- sqrt(v.nri.ev + v.nri.ne)
        #z.nri <- nri/se.nri
        idi <- mean(d[a]) - mean(d[b])
        var.ev <- var(d[a])/na
        var.ne <- var(d[b])/nb
        se.idi <- sqrt(var.ev + var.ne)
        
        nri <- c(nri_cat, se.nri_cat,
                 nri.ev_cat, se.nri.ev_cat,
                 nri.ne_cat, se.nri.ne_cat,
                 nri, se.nri,
                 idi, se.idi)
      })
    reclass_table <-
      do.call("rbind", reclass)
    
    NRI_cat <-
      RR_diff_prop(reclass_table[, 1],
                   reclass_table[, 2])
    NRI_ev <-
      RR_diff_prop(reclass_table[, 3],
                   reclass_table[, 4])
    NRI_ne <-
      RR_diff_prop(reclass_table[, 5],
                   reclass_table[, 6])
    NRI_cont <-
      RR_diff_prop(reclass_table[, 7],
                   reclass_table[, 8])
    IDI <-
      RR_diff_prop(reclass_table[, 9],
                   reclass_table[, 10])
    nri <-
      round(data.frame(rbind(NRI_cat,
                             NRI_ev,
                             NRI_ne,
                             NRI_cont,
                             IDI)), 5)
    names(nri) <-
      c("NRI", "SE", "95% CI LO", "95% CI UP")
    return(nri)
}