\name{HPDI}
\alias{PCI}
\alias{HPDI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Confidence/credible intervals from samples}
\description{
  These functions compute highest posterior density (HPDI) and percentile (PCI) confidence intervals, using samples from the posterior density.
}
\usage{
HPDI( samples , prob=0.95 )
PCI( samples , prob=0.95 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samples}{Vector of parameter values}
  \item{prob}{interval probability mass}
}
\details{
  Highest Posterior Density Intervals (HPDI) are calculated by \code{\link{HPDinterval}} in the \code{coda} package.
  
  Percentile intervals (PCI) use \code{\link{quantile}} and assign equal mass to each tail.
}
\value{
}
\references{McElreath 2011, Statistical Rethinking.}
\author{Richard McElreath}
\seealso{\code{\link{HPDinterval}}}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

