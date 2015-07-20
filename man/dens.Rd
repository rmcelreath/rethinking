\name{dens}
\alias{dens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Density plots}
\description{
  Convenient interface for plotting density estimates.
}
\usage{
dens( x , adj=0.5 , norm.comp=FALSE , main="" , show.HPDI=FALSE , show.zero=FALSE , rm.na=TRUE , add=FALSE , ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector of values to construct density from, or data frame. If \code{x} is a data frame, then \code{dens} plots a grid of densities, one for each column in \code{x}.}
  \item{adj}{width of density kernal.}
  \item{norm.comp}{If \code{TRUE}, overlays normal density comparison.}
  \item{show.HPDI}{If a numeric value, displays HPDI of same width. For example, \code{show.HPDI=0.95} shows a 95 percent HPDI inside the density.}
  \item{show.zero}{If \code{TRUE}, draws a vertical line at location of zero on horizonal axis.}
  \item{rm.na}{If \code{TRUE}, removes \code{NA}s before computing density}
  \item{add}{If \code{TRUE}, overlays density on an existing plot}
  \item{...}{Other parameters to pass to \code{density}, which constructs the density estimates.}
}
\details{
  This function merely provides a convenient interface for plotting density estimates produced by \code{density}. It handles both single vectors and multiple vectors, contained within a data frame.
  
  Highest Posterior Density Intervals (HPDI) are calculated by \code{HPDinterval} in the \code{coda} package.
}
\value{
}
\references{}
\author{Richard McElreath}
\seealso{\code{\link{density}}, \code{\link{HPDinterval}}}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

