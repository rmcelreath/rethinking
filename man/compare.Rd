\name{compare}
\alias{compare}\alias{ICweights}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compare fit models using WAIC or DIC}
\description{
  Returns a table of model comparison statistics, by default focused on WAIC.
}
\usage{
compare( ... , n=1e3 , sort="WAIC" , func=WAIC , WAIC=TRUE , refresh=0 , log_lik="log_lik" )
ICweights( dev )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{A series of fit models, separated by commas}
  \item{n}{Number of samples from posterior to use in computing WAIC/DIC}
  \item{sort}{Sort table by ascending values in named column}
  \item{func}{Function to use in computing criteria for comparison}
  \item{WAIC}{Deprecated: If \code{TRUE}, uses \code{func} for comparison. If \code{FALSE}, uses DIC.}
  \item{refresh}{Progress display update interval. 0 suppresses display.}
  \item{dev}{Vector of values to use in computing model weights}
  \item{log_lik}{Name of log likelihood list in posterior}
}
\details{
  This function computes WAIC and optionally DIC values for fit models and returns a table sorted by ascending values. Each row in this table is a model, and the various columns provide WAIC, effective numbers of parameters, model weights, and standard errors.
  
  A \code{plot} method is supported, for graphic display of the information criteria.
}
\value{
An object of class \code{compareIC} with slots \code{output} (table of results) and \code{dSE} (matrix of standard errors of differences in IC between pairs of models).
}
\references{}
\author{Richard McElreath}
\seealso{}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

