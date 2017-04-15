
#' @export
fitted.map <- function(object, n = 1L, data = NULL, post = NULL, ...,
                       na.rm = TRUE) {
  if (is.null(data)) data <- object@data

  if (is.null(post)) {
    if (n == 1) {
      apply(link(object, post = coef(object), data = data, ...),
            2, mean, na.rm = na.rm)
    } else {
      apply(link(object, n = n, data = data, ...), 2, mean, na.rm = na.rm)
    }
  } else {
    apply(link(object, post = post, data = data, ...), 2, mean, na.rm = na.rm)
  }
}

#' @export
residuals.map <-
  function(object, n = 1L, data = NULL, fits = NULL, post = NULL, ..., na.rm = TRUE,
           silent = !verbose, verbose = TRUE){

  if (is.null(data)) {
    data <- object@data
  }

  if (is.null(fits)) {
    fits <- fitted(object, n=n, data = data, post = post, ..., na.rm = na.rm)
  }
  outcome_name <- object@formula[[1]][[2]]
  if (!silent) {
    message("Residuals computed assuming outcome = ", as.character(outcome_name))
  }
  eval(outcome_name, data) - fits
}