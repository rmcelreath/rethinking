#' Class \code{map2stan}: fitted map2stan Stan model
#'
#' This object contains an object of class \code{stanfit}, return by
#' \code{\link{stan}}, along with other information produced by a
#' \code{\link{map2stan}} model fit.
#'
#'
#' @name map2stan-class
#' @aliases map2stan-class show,map2stan-method pairs,map2stan-method
#' plot,map2stan-method extract.samples,map2stan-method
#' stancode,map2stan-method WAIC,map2stan-method DIC,map2stan-method
#' deviance,map2stan-method sim,map2stan-method logLik,map2stan-method
#' nobs,map2stan-method summary,map2stan-method vcov,map2stan-method
#' @docType class
#' @section Slots: \describe{ \item{list("call")}{Original function call that
#' produced the fit}\item{:}{Original function call that produced the fit}
#' \item{list("model")}{Text of the Stan model used in sampling}\item{:}{Text
#' of the Stan model used in sampling} \item{list("stanfit")}{Object of class
#' \code{stanfit}}\item{:}{Object of class \code{stanfit}}
#' \item{list("coef")}{Vector of posterior means}\item{:}{Vector of posterior
#' means} \item{list("vcov")}{Diagonal matrix with parameter
#' variances}\item{:}{Diagonal matrix with parameter variances}
#' \item{list("data")}{Data used in sampling}\item{:}{Data used in sampling}
#' \item{list("start")}{Initial values used in sampling}\item{:}{Initial values
#' used in sampling} \item{list("pars")}{Vector of parameters returned by
#' sampling}\item{:}{Vector of parameters returned by sampling}
#' \item{list("formula")}{Original \code{alist} of formulas used to define
#' model}\item{:}{Original \code{alist} of formulas used to define model}
#' \item{list("formula_parsed")}{Parsed formula list. Used during compilation
#' and useful for debugging and class methods.}\item{:}{Parsed formula list.
#' Used during compilation and useful for debugging and class methods.} }
#' @seealso \code{\link{map2stan}}
#' @keywords classes
#'
#'
NULL



