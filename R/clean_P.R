#' Function to clean variables
#' 
#' @author Martijn Heymans, 2020
#' @keywords internal
#'
#' @export
clean_P <-
  function(variable){
    variable <-
      gsub(paste(c("factor", "[()]", 
                   "rcs", ",", " ", 
                   c(3:7)), collapse = "|"),
           "", variable)
    return(variable)
}