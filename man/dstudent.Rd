\name{dstudent}
\alias{dstudent}\alias{rstudent}
\title{Non-standard Student's t distribution}
\description{
  Functions for computing density and producing random samples from a non-standardized Student's t distribution.
}
\usage{
dstudent( x, nu = 2, mu = 0, sigma = 1, log = FALSE )
rstudent( n, nu = 2, mu = 0, sigma = 1 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Values to compute densities of}
  \item{nu}{Degrees of freedom (tail shape)}
  \item{mu}{Location of distribution (mean)}
  \item{sigma}{Scale of distribution}
  \item{log}{If \code{TRUE}, returns log-density instead of density}
  \item{n}{Number of random observations to sample}
}
\details{
  These functions provide density and random number calculations for Student's t distribution, translated and scaled by mean \code{mu} and scale \code{sigma}. Note that \code{sigma} is not the distribution's standard deviation, unless \code{nu} is very large.
}
\references{}
\author{Richard McElreath}
\seealso{}
\examples{
\dontrun{
library(rethinking)
data(WaffleDivorce)
d <- WaffleDivorce
d$A <- scale( d$MedianAgeMarriage )
d$D <- scale( d$Divorce )
d$M <- scale( d$Marriage )

m5.3b <- quap(
    alist(
        D ~ dstudent( 2 , mu , sigma ) ,
        mu <- a + bM*M + bA*A ,
        a ~ dnorm( 0 , 0.2 ) ,
        bM ~ dnorm( 0 , 0.5 ) ,
        bA ~ dnorm( 0 , 0.5 ) ,
        sigma ~ dexp( 1 )
    ) , data = d )
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

