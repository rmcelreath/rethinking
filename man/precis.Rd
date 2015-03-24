\name{precis}
\alias{precis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Precis of model fit}
\description{
  Displays concise parameter estimate information for an existing model fit.
}
\usage{
precis( model , depth=1 , pars , ci=TRUE , prob=0.89 , 
    corr=FALSE , digits=2 , warn=TRUE )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{Fit model object}
  \item{depth}{If \code{1}, suppresses vectors and matrices of parameters. If \code{2}, displays all parameters}
  \item{pars}{Optional character vector of parameter names to display}
  \item{ci}{Show quadratic estimate confidence intervals}
  \item{prob}{Width of posterior intervals}
  \item{corr}{If \code{TRUE}, show correlations among parameters in output}
  \item{digits}{Number of decimal places to display in output}
  \item{warn}{If \code{TRUE}, warns about various things}
}
\details{
  Creates a table of estimates and standard errors, with optional confidence intervals and parameter correlations. Posterior intervals are quadratic estimates, derived from standard deviations, unless the model uses samples from the posterior distribution, in which case \code{\link{HPDI}} is used instead.
  
  Can also provide expected value, standard deviation, and HPDI columns for a data frame.
}
\value{
  A data frame with a row for each parameter.
}
\references{}
\author{Richard McElreath}
\seealso{}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

