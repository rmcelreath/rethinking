\name{dzipois}
\alias{dzipois}\alias{rzipois}
\title{Zero-inflated Poisson probability density}
\description{
  Functions for computing density and producing random samples from a zero-inflated Poisson probability distribution.
}
\usage{
dzipois( x , p , lambda , log=FALSE )
rzipois( n , p , lambda )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Integer values to compute densities or probabilies of}
  \item{p}{Probability of a zero from bernoulli process}
  \item{lambda}{Poisson rate parameter (mean)}
  \item{log}{If \code{TRUE}, returns log-probability instead of probability}
}
\details{
  These functions provide density and random number calculations for zero-inflated Poisson observations. This distribution is defined as a finite mixture of zeros and Poisson values, where \code{p} determines the weight of the pure zeros. As such, the probability of a zero is \code{p + (1-p)exp(-lambda)}. And the probability of a non-zero value \code{x} is \code{(1-p)lambda^x exp(-lambda)/x!}.
  
  \code{dzipois} does its calculations in a way meant to reduce numerical issues with probabilities close to 0 and 1. As a result, it's not very fast. 
}
\references{}
\author{Richard McElreath}
\seealso{}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

