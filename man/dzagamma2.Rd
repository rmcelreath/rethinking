\name{dzagamma2}
\alias{dzagamma2}
\title{Zero-augmented gamma probability density}
\description{
  Function for computing density from a zero-augmented gamma probability distribution.
}
\usage{
dzagamma2( x , prob , mu , scale , log=FALSE )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Values to compute densities for}
  \item{prob}{Probability of a zero}
  \item{mu}{Mean of gamma distribution}
  \item{scale}{Scale parameter (same as \code{scale} in \code{\link{dgamma}})}
  \item{log}{If \code{TRUE}, returns log-density instead of density}
}
\details{
  This distribution is defined as a finite mixture of zeros and strictly positive gamma distributed values, where \code{prob} determines the weight of the zeros. As such, the probability of a zero is \code{prob}, and the probability of a non-zero value \code{x} is \code{(1-prob)*dgamma( x , mu/scale , scale )}.
}
\references{}
\author{Richard McElreath}
\seealso{}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

