\name{dzibinom}
\alias{dzibinom}\alias{rzibinom}
\title{Zero-inflated binomial probability density}
\description{
  Functions for computing density and producing random samples from a zero-inflated binomial probability distribution.
}
\usage{
dzibinom( x , p_zero , size , prob , log=FALSE )
rzibinom( n , p_zero , size , prob )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Integer values to compute densities or probabilies of}
  \item{p_zero}{Probability of a zero from bernoulli process}
  \item{size}{Number of binomial trials}
  \item{prob}{Probability of success in binomial trial}
  \item{log}{If \code{TRUE}, returns log-probability instead of probability}
}
\details{
  These functions provide density and random number calculations for zero-inflated binomial observations. This distribution is defined as a finite mixture of zeros and binomial values, where \code{p_zero} determines the weight of the pure zeros. As such, the probability of a zero is \code{p_zero + (1-p_zero)(1-prob)^size}. And the probability of a non-zero value \code{x} is \code{(1-p_zero) choose(size,x) prob^x (1-prob)^(size-x)}.
  
  \code{dzibinom} does its calculations in a way meant to reduce numerical issues with probabilities close to 0 and 1. As a result, it's not very fast. 
}
\references{}
\author{Richard McElreath}
\seealso{}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

