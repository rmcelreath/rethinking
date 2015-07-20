\name{sample.qa.posterior}
\alias{sample.naive.posterior}
\alias{sample.qa.posterior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Samples from quadratic posterior densities of models}
\description{
  Samples from the posterior density of a fit model or models, assuming multivariate normal density. 
}
\usage{
sample.qa.posterior(model, n = 10000, clean.names = TRUE, model.weights =
                 "AICc", nobs = 0, add.names = FALSE, fill.na = 0,
                 verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{A fit model object}
  \item{models}{A list of fit models of the same class}
  \item{n}{Number of samples to draw from joint posterior}
  \item{model.weights}{If passing a list of models, method for computing posterior probability of each model family. Can be "AIC","AICc","BIC" or a vector of numeric weights.}
  \item{nobs}{Number of observations used to fit model or all models in list. Sometimes needed for \code{model.weights} values, like \code{AICc}.}
  \item{add.names}{Adds a column of model names, when passing a list of models}
  \item{fill.na}{Fills missing values with 0, by default, for model families that do not contain a given parameter. Useful for linear models. Hazardous for non-linear ones.}
  \item{verbose}{If \code{TRUE}, prints various debugging information}
}
\details{
  This is a legacy function and is no longer supported nor unit tested.

  This function provides a way to draw parameter values from a multivariate normal posterior density, estimated from the maximum a posterieri (MAP) estimates and variance-covariance (\code{vcov}) of a fit model or models.
  
  When passing a single fit model object, the function returns a data frame in which each row is a sample and each column is a parameter.
  
  When passing a list of fit model objects, the function returns a data frame containing samples from the joint posterior across model families. The fraction of rows drawn from a specific model family is determined by the \code{model.weights} parameter. BIC, AIC, or AICc are used to compute approximate predictive probabilities of each model family, and the total samples \code{n} is proportioned according to these estimates. The user can also supply a numeric vector of model weights, computed by any method. This vector should sum to 1.
}
\value{
}
\references{}
\author{Richard McElreath}
\seealso{\code{\link{mvrnorm}}}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

