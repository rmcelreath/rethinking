\name{dgampois}
\alias{dgampois}\alias{rgampois}
\title{Gamma-Poisson probability density}
\description{
  Functions for computing density and producing random samples from a gamma-Poisson (negative-binomial) probability distribution.
}
\usage{
dgampois( x , mu , scale , log=FALSE )
rgampois( n , mu , scale )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Integer values to compute probabilies of}
  \item{mu}{Mean of gamma distribution}
  \item{scale}{Scale parameter of gamma distribution}
  \item{log}{If \code{TRUE}, returns log-probability instead of probability}
  \item{n}{Number of random observations to sample}
}
\details{
  These functions provide density and random number calculations for gamma-Poisson observations. These functions use \code{dnbinom} and \code{rnbinom} internally, but convert the parameters from the \code{mu} and \code{scale} form. The \code{size} parameter of the negative-binomial is defined by \code{mu/scale}, and the \code{prob} parameter of the negative-binomial is the same as \code{1/(1+scale)}.
}
\references{}
\author{Richard McElreath}
\seealso{}
\examples{
\dontrun{
data(Hurricanes)

# map model fit
# note exp(log_scale) to constrain scale to positive reals
m <- map(
    alist(
        deaths ~ dgampois( mu , exp(log_scale) ),
        log(mu) <- a + b*femininity,
        a ~ dnorm(0,100),
        b ~ dnorm(0,1),
        log_scale ~ dnorm(1,10)
    ),
    data=Hurricanes )

# map2stan model fit
# constraint on scale is passed via contraints list
m.stan <- map2stan(
    alist(
        deaths ~ dgampois( mu , scale ),
        log(mu) <- a + b*femininity,
        a ~ dnorm(0,100),
        b ~ dnorm(0,1),
        scale ~ dcauchy(0,2)
    ),
    data=Hurricanes,
    constraints=list(scale="lower=0"),
    start=list(scale=2) )
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

