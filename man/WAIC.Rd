\name{WAIC}
\alias{WAIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Widely Applicable Information Criterion}
\description{
  Computes WAIC for a \code{map2stan} model fit.
}
\usage{
WAIC( object , n=1000 , refresh=0.1 , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class \code{map2stan}}
  \item{n}{Number of samples to use in computing WAIC. Set to \code{n=0} to use all samples in \code{map2stan} fit}
  \item{refresh}{Refresh interval for progress display. Set to \code{refresh=0} to suppress display.}
  \item{...}{Other parameters to pass to someone}
}
\details{
  This function uses the samples and model definition from a \code{map2stan} fit to compute the Widely Applicable Information Criterion, WAIC. WAIC is an estimate of out-of-sample relative K-L divergence, and it is defined as:

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
\seealso{\code{\link{map2stan}}, \code{\link{link}}}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

