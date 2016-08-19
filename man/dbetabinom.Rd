\name{dbetabinom}
\alias{dbetabinom}\alias{rbetabinom}
\title{Beta-binomial probability density}
\description{
  Functions for computing density and producing random samples from a beta-binomial probability distribution.
}
\usage{
dbetabinom( x , size , prob , theta , shape1 , shape2 , log=FALSE )
rbetabinom( n , size , prob , theta , shape1 , shape2 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Integer values to compute probabilies of}
  \item{size}{Number of trials}
  \item{prob}{Average probability of beta distribution}
  \item{theta}{Dispersion of beta distribution}
  \item{shape1}{First shape parameter of beta distribution (alpha)}
  \item{shape2}{Second shape parameter of beta distribution (beta)}
  \item{log}{If \code{TRUE}, returns log-probability instead of probability}
  \item{n}{Number of random observations to sample}
}
\details{
  These functions provide density and random number calculations for beta-binomial observations. The \code{dbetabinom} code is based on Ben Bolker's original code in the \code{emdbook} package. 
  
  Either \code{prob} and \code{theta} OR \code{shape1} and \code{shape2} must be provided. The two parameterizations are related by shape1 = prob * theta, shape2 = (1-prob) * theta.
  
  The \code{rbetabinom} function generates random beta-binomial observations by using both \code{\link{rbeta}} and \code{\link{rbinom}}. It draws \code{n} values from a beta distribution. Then for each, it generates a random binomial observation.
}
\references{}
\author{Richard McElreath}
\seealso{}
\examples{
\dontrun{
data(reedfrogs)
reedfrogs$pred_yes <- ifelse( reedfrogs$pred=="pred" , 1 , 0 )

# map model fit
# note exp(log_theta) to constrain theta to positive reals
m <- map(
    alist(
        surv ~ dbetabinom( density , p , exp(log_theta) ),
        logit(p) <- a + b*pred_yes,
        a ~ dnorm(0,10),
        b ~ dnorm(0,1),
        log_theta ~ dnorm(1,10)
    ),
    data=reedfrogs )

# map2stan model fit
# constraint on theta is passed via contraints list
m.stan <- map2stan(
    alist(
        surv ~ dbetabinom( density , p , theta ),
        logit(p) <- a + b*pred_yes,
        a ~ dnorm(0,10),
        b ~ dnorm(0,1),
        theta ~ dcauchy(0,1)
    ),
    data=reedfrogs,
    constraints=list(theta="lower=0") )
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

