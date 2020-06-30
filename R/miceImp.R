#' Wrapper function around mice
#' 
#'@author Martijn Heymans, 2020
#'@keywords internal
#'
#'@export
miceImp <- function(data, ...) {
  Imp <- mice::mice(data, ...)
}

