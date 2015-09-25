\name{WAIC}
\alias{WAIC}\alias{DIC}\alias{LOO}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Information Criteria and Pareto-Smoothed Importance Sampling Cross-Validation}
\description{
  Computes WAIC, DIC, and PSIS-LOO for \code{map} and \code{map2stan} model fits. In addition, WAIC and PSIS-LOO can be calculated for Stan model fits (see details).
}
\usage{
WAIC( object , n=1000 , refresh=0.1 , pointwise=FALSE , ... )
LOO( object , n=1000 , refresh=0.1 , pointwise=FALSE , ... )
DIC( object , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class \code{map} or \code{map2stan}}
  \item{n}{Number of samples to use in computing WAIC. Set to \code{n=0} to use all samples in \code{map2stan} fit}
  \item{refresh}{Refresh interval for progress display. Set to \code{refresh=0} to suppress display.}
  \item{pointwise}{If \code{TRUE}, return a vector of WAIC values for each observation. Useful for computing standard errors.}
  \item{...}{Other parameters to pass to specific methods}
}
\details{
  These functions use the samples and model definition from a \code{map} or \code{map2stan} fit to compute the Widely Applicable Information Criterion (WAIC), Deviance Information Criterion (DIC), or Pareto-smoothed importance-sampling leave-one-out cross-validation estimate (PSIS-LOO). 
  
  WAIC is an estimate of out-of-sample relative K-L divergence (KLD), and it is defined as:

  \deqn{WAIC = -2(lppd - pWAIC)}

  Components \code{lppd} (log pointwise predictive density) and \code{pWAIC} (the effective number of parameters) are reported as attributes. See Gelman et al 2013 for definitions and formulas. This function uses the variance definition for \code{pWAIC}.

  PSIS-LOO is another estimate of out-of-sample relative K-L divergence. It is computed by the \code{loo} package. See Vehtari et al 2015 for definitions and computation.

  In practice, WAIC and PSIS-LOO are extremely similar estimates of KLD.

  Both WAIC and LOO have methods for \code{stanfit} models, provided the posterior contains a log-likelihood matrix (samples on rows, observations on columns) named \code{log_lik}. See example.
}
\value{
}
\references{
Watanabe, S. 2010. Asymptotic equivalence of Bayes cross validation and Widely Applicable Information Criterion in singular learning theory. Journal of Machine Learning Research 11:3571-3594.

Gelman, A., J. Hwang, and A. Vehtari. 2013. Understanding predictive information criteria for Bayesian models.

Vehtari, A., A. Gelman, and J. Gabry. 2015. Efficient implementation of leave-one-out cross-validation and WAIC for evaluating fitted Bayesian models.
}
\author{Richard McElreath}
\seealso{\code{\link{map}}, \code{\link{map2stan}}, \code{\link{link}}, \code{\link{loo}}}
\examples{
\dontrun{
library(rethinking)
data(chimpanzees)
d <- chimpanzees
dat <- list(
    y = d$pulled_left,
    prosoc = d$prosoc_left,
    condition = d$condition,
    N = nrow(d)
)

m1s_code <- '
data{
    int<lower=1> N;
    int y[N];
    int prosoc[N];
}
parameters{
    real a;
    real bP;
}
model{
    vector[N] p;
    bP ~ normal( 0 , 1 );
    a ~ normal( 0 , 10 );
    for ( i in 1:N ) {
        p[i] <- a + bP * prosoc[i];
    }
    y ~ binomial_logit( 1 , p );
}
generated quantities{
    vector[N] p;
    vector[N] log_lik;
    for ( i in 1:N ) {
        p[i] <- a + bP * prosoc[i];
        log_lik[i] <- binomial_logit_log( y[i] , 1 , p[i] );
    }
}
'

m1s <- stan( model_code=m1s_code , data=dat , chains=2 , iter=2000 )

WAIC(m1s)

LOO(m1s)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

