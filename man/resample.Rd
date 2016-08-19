\name{resample}
\alias{resample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Resample map2stan fit}
\description{
  Sample from a new chain or chains, using a previous \code{map2stan} fit object.
}
\usage{
resample( object , iter=1e4 , warmup=1000 , chains=1 , cores=1 , 
    DIC=TRUE , WAIC=TRUE , rng_seed , data , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class \code{map2stan}}
  \item{iter}{Number of sampling iterations, including warmup}
  \item{warmup}{Number of adaptation steps}
  \item{chains}{Number of independent chains}
  \item{cores}{Number of cores to distribute chains across}
  \item{DIC}{If \code{TRUE}, computes DIC after sampling}
  \item{WAIC}{If \code{TRUE}, computes WAIC after sampling}
  \item{rng_seed}{Optional seed to use for all chains. When missing, a random seed is chosen and used for all chains.}
  \item{...}{Other parameters to pass to \code{stan}}
}
\details{
  This function is a convenience for drawing more samples from an initial \code{map2stan} fit.
  
  When \code{cores} is set greater than 1, either \code{\link{mclapply}} (on a unix system) or \code{\link{parLapply}} (on a Windows system) is used to run the chains, distributing them across processor cores. The results are automatically recombined with \code{\link{sflist2stanfit}}. 
}
\value{
  An object of class \code{map2stan}, holding the new samples, as well as all of the original formulas and data for the model.
}
\references{}
\author{Richard McElreath}
\seealso{\code{\link{map2stan}}, \code{\link{mclapply}}, \code{\link{sflist2stanfit}}}
\examples{
\dontrun{
data(Trolley)
d <- Trolley
d2 <- list(
    y=d$response,
    xA=d$action,
    xI=d$intention,
    xC=d$contact,
    id=as.integer(d$id)
)
Nid <- length(unique(d2$id))

# ordered logit regression with varying intercepts
m.init <- map2stan(
    alist(
        y ~ dordlogit( phi , cutpoints ),
        phi <- aj + bA*xA + bI*xI + bC*xC,
        c(bA,bI,bC) ~ dnorm(0,1),
        aj[id] ~ dnorm(0,sigma_id),
        sigma_id ~ dcauchy(0,2.5),
        cutpoints ~ dcauchy(0,2.5)
    ),
    data=d2 ,
    start=list(
        bA=0,bI=0,bC=0,
        cutpoints=c(-2,-1.7,-1,-0.2,0.5,1.3),
        aj=rep(0,Nid),sigma_id=1
    ),
    types=list(cutpoints="ordered") , 
    iter=2
)

# Note: parallel chains won't work on Windows
m <- resample( m.init , chains=3 , cores=3 , warmup=1000 , iter=3000 )
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

