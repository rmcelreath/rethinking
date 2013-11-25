\name{compare}
\alias{compare}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compare fit models using AICc/DIC/BIC}
\description{
  Returns a table of model comparison statistics
}
\usage{
compare( ... , nobs=NULL , sort="AICc" , BIC=FALSE , DIC=FALSE , delta=TRUE , 
	DICsamples=1e4 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{A series of fit models, separated by commas}
  \item{nobs}{Number of observations to use in calculating AICc/BIC}
  \item{sort}{Sort table by ascending values in named column}
  \item{BIC}{Compute BIC as well}
  \item{DIC}{Compute DIC as well. See details.}
  \item{delta}{Include differences between AICc/BIC/DIC values}
  \item{DICsamples}{Number of posterior samples to use in computing DIC}
}
\details{
  This function computes AICc and optionally BIC and DIC values for fit models and returns a table sorted by ascending values. Each row in this table is a model, and the various columns provide AICc, BIC, DIC, numbers of parameters, and model weights.

  At the moment, DIC only works for models fit with \code{map}. DIC is computing using the original \code{Dhat + 2(Dbar-Dhat)} formulation, not the variance formulation. Samples are drawn from the posterior, assuming it is multivariate gaussian and defined by the variance-covariance matrix of the MAP estimates. The column \code{pD} in the results is the effective number of parameters, \code{Dbar-Dhat}.
}
\value{
}
\references{McElreath 2011, Statistical Rethinking.}
\author{Richard McElreath}
\seealso{}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

