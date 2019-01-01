\name{sim_train_test}
\alias{sim_train_test}
\alias{sim.train.test}
\title{Simulate in-sample and out-of-sample model performance}
\description{
    Simulates in-sample and out-of-sample model performance for simple \code{quap} models, returning lppd in and out of sample, WAIC, LOOIC, and LOOCV.
}
\usage{
sim_train_test(N = 20, k = 3, rho = c(0.15, -0.4), b_sigma = 100, 
    WAIC = FALSE, LOOCV=FALSE , LOOIC=FALSE , cv.cores=1 , 
    return_model=FALSE )

# old function from 1st edition
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
  \item{LOOIC}{If \code{TRUE}, returns LOOIC as produced by \code{LOO}}
  \item{LOOCV}{If \code{TRUE}, returns LOOCV as produced by \code{cv_quap}}
  \item{devbar}{If \code{TRUE}, returns the average deviance in-sample}
  \item{devbarout}{If \code{TRUE}, returns average deviance out-of-sample}
  \item{cv.cores}{Number or cores to use for cross-validation}
  \item{return_model}{If \code{TRUE}, includes fit model in result}
}
\details{
  This function simulates Gaussian data and then fits linear regression models to it, returning the lppd of the fit as produced by \code{lppd} (training, in-sample) and the deviance on a new sample, computed using the posterior from the training sample. 
}
\value{
}
\references{
}
\author{Richard McElreath}
\seealso{ \code{\link{lppd}}, \code{\link{WAIC}}, \code{\link{LOO}}, \code{\link{cv_quap}} }
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

