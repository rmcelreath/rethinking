\name{pairs.map2stan}
\alias{pairs.map2stan}
\alias{pairs.map}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Trace plot of map2stan chains}
\description{
  Shows \code{pairs} plots for Stan samples produced by a \code{map} or \code{map2stan} fit. 
}
\usage{
pairs(x, n=500 , alpha=0.7 , cex=0.7 , pch=16 , adj=1 , pars , ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{map or map2stan model fit}
  \item{n}{Number of samples to show in scatter plots}
  \item{alpha}{alpha transparency of scatter plots}
  \item{cex}{Character expansion factor for scatter plots}
  \item{pch}{Point type for scatter plots}
  \item{adj}{\code{adj} argument for \code{\link{density}} estimate}
  \item{pars}{Character vector of parameters to display}
  \item{...}{Additional plot arguments}
}
\details{
  This is the default pairs plot method for both \code{\link{map}} and \code{\link{map2stan}} model fits. 
}
\value{
}
\references{}
\author{Richard McElreath}
\seealso{}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

