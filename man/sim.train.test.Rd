\name{sim.train.test}
\alias{sim.train.test}
\title{Simulate in-sample and out-of-sample deviance}
\description{
  Simulates posterior observations for \code{map} and \code{map2stan} model fits.
}
\usage{
sim.train.test( N=20 , k=3 , rho=c(0.15,-0.4) , b_sigma=100 , 
    DIC=FALSE , WAIC=FALSE, devbar=FALSE , devbarout=FALSE )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{Number of cases in simulated data}
  \item{k}{Number of parameters in model to fit to data}
  \item{rho}{Vector of correlations between predictors and outcome, in simulated data}
  \item{b_sigma}{Standard deviation of beta-coefficient priors}
  \item{DIC}{If \code{TRUE}, returns DIC}
  \item{WAIC}{If \code{TRUE}, returns WAIC}
  \item{devbar}{If \code{TRUE}, returns the average deviance in-sample}
  \item{devbarout}{If \code{TRUE}, returns average deviance out-of-sample}
}
\details{
  This function simulates Gaussian data and then fits linear regression models to it, returning the deviance of the fit (training, in-sample) and the deviance on a new sample, computed using the MAP estimates from the training sample. 
}
\value{
}
\references{
}
\author{Richard McElreath}
\seealso{}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

