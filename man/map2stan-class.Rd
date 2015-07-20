\name{map2stan-class} 
\Rdversion{1.1}

\docType{class}
\alias{map2stan-class}
\alias{show,map2stan-method}
\alias{pairs,map2stan-method}
\alias{plot,map2stan-method}
\alias{extract.samples,map2stan-method}
\alias{stancode,map2stan-method}
\alias{WAIC,map2stan-method}
\alias{DIC,map2stan-method}
\alias{deviance,map2stan-method}
\alias{sim,map2stan-method}
\alias{logLik,map2stan-method}
\alias{nobs,map2stan-method}
\alias{summary,map2stan-method}
\alias{vcov,map2stan-method}

\title{Class \code{map2stan}: fitted map2stan Stan model}

\description{
  This object contains an object of class \code{stanfit}, return by \code{\link{stan}}, along with other information produced by a \code{\link{map2stan}} model fit.
} 
\section{Slots}{
  \describe{
    \item{\code{call}:}{Original function call that produced the fit}
    \item{\code{model}:}{Text of the Stan model used in sampling}
    \item{\code{stanfit}:}{Object of class \code{stanfit}}
    \item{\code{coef}:}{Vector of posterior means}
    \item{\code{vcov}:}{Diagonal matrix with parameter variances}
    \item{\code{data}:}{Data used in sampling}
    \item{\code{start}:}{Initial values used in sampling}
    \item{\code{pars}:}{Vector of parameters returned by sampling}
    \item{\code{formula}:}{Original \code{alist} of formulas used to define model}
    \item{\code{formula_parsed}:}{Parsed formula list. Used during compilation and useful for debugging and class methods.}
  }
}
\section{Methods}{
  \describe{
    \item{\code{show}}{\code{signature(object = "map2stan")}: print the default summary for the model.}
    \item{\code{stancode}}{Displays the \code{model} slot.}
    \item{\code{WAIC}}{Computes WAIC}
  }
}
\references{
}
\seealso{
  \code{\link{map2stan}} 
}
\examples{
}

\keyword{classes}