\name{HPDI}
\alias{PI}
\alias{HPDI}
\alias{PCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Confidence/credible intervals from samples}
\description{
  These functions compute highest posterior density (HPDI) and percentile (PI) intervals, using samples from a posterior density or simulated outcomes.
}
\usage{
HPDI( samples , prob=0.89 )
PI( samples , prob=0.89 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samples}{Vector of parameter values}
  \item{prob}{interval probability mass}
}
\details{
  Highest Posterior Density Intervals (HPDI) are calculated by \code{\link{HPDinterval}} in the \code{coda} package.
  
  Percentile intervals (PI) use \code{\link{quantile}} and assign equal mass to each tail.
}
\value{
}
\references{}
\author{Richard McElreath}
\seealso{\code{\link{HPDinterval}}}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

