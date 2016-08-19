\name{extract.samples}
\alias{extract.samples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Collect posterior samples from a fit model}
\description{
  Extracts or draws posterior samples for fit models. 
}
\usage{
extract.samples( object , n=10000 , clean.names=TRUE , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Fit model to extract samples from}
  \item{n}{Number of samples to simulate}
  \item{clean.names}{When \code{TRUE}, removes parentheses from parameters named \code{(Intercept)}.}
  \item{...}{Other parameters to pass to descendent functions (when defined)}
}
\details{
  For \code{map2stan} models, this function returns cleaned samples already contained in the object. These samples are cleaned of dimension attributes and the \code{lp__} and \code{dev} traces that are used internally.
  
  For \code{map} and other types, it uses the variance-covariance matrix and coefficients to define a multivariate Gaussian posterior to draw \code{n} samples from.
}
\value{
A named \code{list} (for \code{map2stan}) or \code{data.frame} containing samples for each parameter in the posterior distribution.
}
\references{}
\author{Richard McElreath}
\seealso{\code{\link{mvrnorm}}}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

