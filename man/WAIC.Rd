\name{WAIC}
\alias{WAIC}\alias{DIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Widely Applicable Information Criterion / Deviance Information Criterion}
\description{
  Computes WAIC or DIC for \code{map} and \code{map2stan} model fits.
}
\usage{
WAIC( object , n=1000 , refresh=0.1 , pointwise=FALSE , ... )
DIC( object , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class \code{map} or \code{map2stan}}
  \item{n}{Number of samples to use in computing WAIC. Set to \code{n=0} to use all samples in \code{map2stan} fit}
  \item{refresh}{Refresh interval for progress display. Set to \code{refresh=0} to suppress display.}
  \item{pointwise}{If \code{TRUE}, return a vector of WAIC values for each observation. Useful for computing standard errors.}
  \item{...}{Other parameters to pass to someone}
}
\details{
  These function use the samples and model definition from a \code{map} or \code{map2stan} fit to compute the Widely Applicable Information Criterion, WAIC, or rather the Deviance Information Criterion, DIC. 
  
  WAIC is an estimate of out-of-sample relative K-L divergence, and it is defined as:

  \deqn{WAIC = -2(lppd - pWAIC)}

  Components \code{lppd} (log pointwise predictive density) and \code{pWAIC} (the effective number of parameters) are reported as attributes. See Gelman et al 2013 for definitions and formulas. This function uses the variance definition for \code{pWAIC}.
  
  The function \code{\link{link}} is used internally the compute the values of any linear models.
}
\value{
}
\references{
Watanabe, S. 2010. Asymptotic equivalence of Bayes cross validation and Widely Applicable Information Criterion in singular learning theory. Journal of Machine Learning Research 11:3571-3594.

Gelman, A., J. Hwang, and A. Vehtari. 2013. Understanding predictive information criteria for Bayesian models.
}
\author{Richard McElreath}
\seealso{\code{\link{map}}, \code{\link{map2stan}}, \code{\link{link}}}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

