\name{extract.samples}
\alias{extract.samples}\alias{extract.prior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Collect posterior or prior samples from a map or map2stan model}
\description{
  Extracts or draw samples from fit models.
}
\usage{
extract.samples( object , n=10000 , pars , ... )
extract.prior( object , n=1000 , pars , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Fit model to extract samples from}
  \item{n}{Number of samples to simulate}
  \item{pars}{Character vector of parameters to return}
  \item{...}{Other parameters to pass to descendent functions (when defined)}
}
\details{
  Use \code{extract.samples} and \code{extract.prior} to return lists of samples for posterior and prior distributions, respectively. Methods for \code{extract.samples} are provided for \code{map}, \code{map2stan}, and \code{stanfit} objects. Methods for \code{extract.prior} are provided for \code{map} and \code{map2stan} objects. 

  For \code{map2stan} and \code{stanfit} models, \code{extract.samples} returns cleaned samples already contained in the object. These samples are cleaned of dimension attributes and the \code{lp__}, \code{dev}, and \code{log_lik} traces that are used internally. For \code{map} and other types, it uses the variance-covariance matrix and coefficients to define a multivariate Gaussian posterior to draw \code{n} samples from.

  \code{extract.prior} must simulate draws using the model definition. It attempts to return samples in the same list structure as posterior samples. This makes prior-posterior comparison easier. See examples below.
}
\value{
A named \code{list} (for \code{map2stan} and \code{stanfit}) or \code{data.frame} containing samples for each parameter in the posterior/prior distribution.
}
\references{}
\author{Richard McElreath}
\seealso{\code{\link{mvrnorm}}}
\examples{
data(chimpanzees)

d <- list( 
    pulled_left = chimpanzees$pulled_left ,
    prosoc_left = chimpanzees$prosoc_left ,
    condition = chimpanzees$condition ,
    actor = as.integer( chimpanzees$actor ) ,
    blockid = as.integer( chimpanzees$block )
)

m <- map(
    alist(
        pulled_left ~ dbinom(1,theta),
        logit(theta) <- a + aj[actor] + bp*prosoc_left + bpc*condition*prosoc_left,
        aj[actor] ~ dnorm( 0 , 1 ),
        a ~ dnorm(0,2),
        bp ~ dnorm(0,1),
        bpc ~ dnorm(0,1)
    ) ,
    data=d )

prior <- extract.prior(m,n=1e4)
post <- extract.samples(m)
ps <- par("bty")
par(bty="n")
plot( precis(prior,2) , col.ci="gray" , xlim=c(-3,3.5) , bty="n" )
plot( precis(post,2) , add=TRUE , pch=16 )
par(bty=ps)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

