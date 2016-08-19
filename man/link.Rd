\name{link-methods}
\docType{methods}
\alias{link}
\alias{link-methods}
\alias{link,map2stan-method}
\alias{link,map-method}
\alias{link,lm-method}

\title{Predictions for map and map2stan models}
\description{
  Computes inverse-link linear model values for \code{map} and \code{map2stan} samples.
}
\usage{
link( fit , data , n=1000 , ... )
\S4method{link}{map2stan}( fit , data , n=1000 , post , refresh=0.1 , 
    replace=list() , flatten=TRUE , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{Object of class \code{map} or \code{map2stan}}
  \item{data}{Optional list of data to compute predictions over. When missing, uses data found inside fit object.}
  \item{n}{Number of samples to use}
  \item{post}{Optional samples from posterior. When missing, \code{link} extracts the samples using \code{\link{extract.samples}}.}
  \item{refresh}{Refresh interval for progress display. Set to \code{refresh=0} to suppress display.}
  \item{replace}{Optional named list of samples to replace inside posterior samples. See examples.}
  \item{flatten}{When \code{TRUE}, removes linear model names from result}
  \item{...}{Other parameters to pass to someone}
}
\details{
  This function computes the value of each linear model at each sample for each case in the data. Inverse link functions are applied, so that for example a logit link linear model produces probabilities, using the logistic transform. 

  This function is used internally by \code{\link{WAIC}}, \code{\link{sim}}, \code{\link{postcheck}}, and \code{\link{ensemble}}.
  
  It is possible to replace components of the posterior distribution with simulated values. The \code{replace} argument should be a named \code{list} with replacement values. This is useful for marginalizing over varying effects. See the examples below for an example in which varying intercepts are marginalized this way.
}
\value{
}
\references{
}
\author{Richard McElreath}
\seealso{\code{\link{map}}, \code{\link{map2stan}}, \code{\link{sim}}, \code{\link{ensemble}}, \code{\link{postcheck}}}
\examples{
\dontrun{
library(rethinking)
data(chimpanzees)
d <- chimpanzees
d$recipient <- NULL     # get rid of NAs

# model 4 from chapter 12 of the book
m12.4 <- map2stan( 
    alist(
        pulled_left ~ dbinom( 1 , p ) ,
        logit(p) <- a + a_actor[actor] + (bp + bpC*condition)*prosoc_left ,
        a_actor[actor] ~ dnorm( 0 , sigma_actor ),
        a ~ dnorm(0,10),
        bp ~ dnorm(0,10),
        bpC ~ dnorm(0,10),
        sigma_actor ~ dcauchy(0,1)
    ) ,
    data=d , warmup=1000 , iter=4000 , chains=4 )

# posterior predictions for a particular actor
chimp <- 2
d.pred <- list(
    prosoc_left = c(0,1,0,1),   # right/left/right/left
    condition = c(0,0,1,1),     # control/control/partner/partner
    actor = rep(chimp,4)
)
link.m12.4 <- link( m12.4 , data=d.pred )
apply( link.m12.4 , 2 , mean )
apply( link.m12.4 , 2 , PI )

# posterior predictions marginal of actor
# here we replace the varying intercepts samples 
#   with simulated values

# replace varying intercept samples with simulations
post <- extract.samples(m12.4)
a_actor_sims <- rnorm(7000,0,post$sigma_actor)
a_actor_sims <- matrix(a_actor_sims,1000,7)

# fire up link
# note use of replace list
link.m12.4 <- link( m12.4 , n=1000 , data=d.pred , 
    replace=list(a_actor=a_actor_sims) )

# summarize
apply( link.m12.4 , 2 , mean )
apply( link.m12.4 , 2 , PI )
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

