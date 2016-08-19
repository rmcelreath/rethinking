\name{sim}
\alias{sim}
\alias{sim-methods}
\alias{sim,map2stan-method}
\alias{sim,map-method}
\alias{sim,lm-method}
\title{Simulate posterior observations}
\description{
  Simulates posterior observations for \code{map} and \code{map2stan} model fits.
}
\usage{
sim( fit , data , n=1000 , post , ll , refresh=0.1 , replace , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{Object of class \code{map} or \code{map2stan}}
  \item{data}{Optional list of data to compute predictions over}
  \item{n}{Number of samples to use}
  \item{post}{Samples from posterior. If missing, \code{sim} samples using \code{n}}
  \item{...}{Other parameters to pass to someone}
}
\details{
  This function uses the model definition from a \code{map} or \code{map2stan} fit to simulate outcomes that average over the posterior distribution. It uses \code{\link{link}} internally to process any linear models. Then the correspond random number function is used, as defined by the model's likelihood.
  
  The \code{rethinking} package defines a generic function \code{sim}, so methods can be defined for other model fit classes. I might eventually build methods for \code{lm} and \code{glm}.
}
\value{
}
\references{
}
\author{Richard McElreath}
\seealso{\code{\link{link}},\code{\link{map}},\code{\link{map2stan}}}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

