\name{resample}
\alias{resample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Resample map2stan fit}
\description{
  Sample from a new chain or chains, using a previous \code{map2stan} fit object.
}
\usage{
resample( object , iter=1e4 , warmup=1000 , chains=1 , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class \code{map2stan}}
  \item{iter}{Number of sampling iterations, including warmup}
  \item{warmup}{Number of adaptation steps}
  \item{chains}{Number of independent chains}
  \item{...}{Other parameters to pass to \code{stan}}
}
\details{
  This function is a convenience for drawing more samples from an initial \code{map2stan} fit.
}
\value{
  An object of class \code{map2stan}, holding the new samples, as well as all of the original formulas and data for the model.
}
\references{McElreath 2011, Statistical Rethinking.}
\author{Richard McElreath}
\seealso{\code{\link{map2stan}}}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

