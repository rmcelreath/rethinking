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
  \item{n}{Number of samples to use in computing WAIC}
  \item{refresh}{Refresh interval for progress display. Set to \code{refresh=0} to suppress display.}
  \item{...}{Other parameters to pass to someone}
}
\details{
  This function uses the samples and model definition from a \code{map2stan} fit to compute the Widely Applicable Information Criterion, WAIC. Components \code{lppd} and \code{pWAIC} (the effective number of parameters) are reported as attributes.
  
  The function \code{\link{link}} is used internally the compute the values of any linear models.
}
\value{
}
\references{McElreath 2011, Statistical Rethinking.}
\author{Richard McElreath}
\seealso{\code{\link{map2stan}}, \code{\link{link}}}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

