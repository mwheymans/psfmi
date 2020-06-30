#' Function to clean variables
#' 
#' @param variable Character vector of variable names
#'
#' @author Martijn Heymans, 2020
#' @keywords internal
#'
#' @export
clean_P <-
  function(variable){
    variable <- 
      stringr::str_remove(variable, "factor")
    variable <- 
      stringr::str_remove_all(variable, "[()]") 
    variable <- 
      stringr::str_remove(variable, "rcs") 
    variable <- 
      stringr::str_remove(variable, ",") 
    variable <- 
      stringr::str_squish(variable)
    variable <- 
      stringr::str_remove_all(variable, paste(c(3:7), collapse = "|")) 
    variable <-
      gsub(" ", "", variable)
    return(variable)
  }
