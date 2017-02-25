fitted.map <- function(object, n = 1e3, data = object@data, ..., na.rm = TRUE) {
  apply(link(object, n = n, data = data, ...), 2, mean, na.rm = na.rm)
}

residuals.map <- function(object, n = 1e3, data = object@data, fits = NULL, ..., na.rm = TRUE){
  if (is.null(fits)) {
    fits <- fitted(object, n=n, data = data, ..., na.rm = na.rm)
  }
  outcome_name <- object@formula[[1]][[2]]
  message("Residuals computed assuming outcome = ", as.character(outcome_name))
  eval(outcome_name, data) - fits
}