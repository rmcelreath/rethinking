\name{axis_unscale}
\alias{axis_unscale}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Draw an axis with units on original scale}
\description{
  When plotting a standardized or rescaled variable, this function draws the axis units on the original scale.
}
\usage{
axis_unscale( side = 1, at, orig, factor, ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{side}{Side for axis. 1 is bottom.}
  \item{at}{Locations of tick marks, in original scale values}
  \item{orig}{The variable on original scale. Use this when variable was standardized.}
  \item{factor}{Factor the original variable was multiplied by to get rescaled variable. Use this when rescaling by a reference value, for example dividing by maximum value.}
}
\details{
  This function draws a plot axis with display units on original scale. The typical situation for using this is when an analysis was performed on a standardized or rescaled variable. Plotting the posterior predictions with units on the transformed scale can make interpretation difficult.

  When the variable was standardized (mean subtracted and divided by standard devation) before analysis, use the \code{orig} argument to point to the variable on the original scale.

  When the variable was rescaled (multiplied by a factor to rescale, without relocating zero) before analysis, use the \code{factor} argument.
}
\value{
}
\references{}
\author{Richard McElreath}
\seealso{}
\examples{
sppnames <- c( "afarensis","africanus","habilis","boisei",
    "rudolfensis","ergaster","sapiens")
brainvolcc <- c( 438 , 452 , 612, 521, 752, 871, 1350 )
masskg <- c( 37.0 , 35.5 , 34.5 , 41.5 , 55.5 , 61.0 , 53.5 )
d <- data.frame( species=sppnames , brain=brainvolcc , mass=masskg )
d$mass_std <- (d$mass - mean(d$mass))/sd(d$mass)
d$brain_std <- d$brain / max(d$brain)

m7.2 <- quap(
  alist(
      brain_std ~ dnorm( mu , exp(log_sigma) ),
      mu <- a + b[1]*mass_std + b[2]*mass_std^2,
      a ~ dnorm( 0.5 , 1 ),
      b ~ dnorm( 0 , 10 ),
      log_sigma ~ dnorm( 0 , 1 )
  ), data=d , start=list(b=rep(0,2)) )

plot( d$brain_std ~ d$mass_std , xaxt="n" , yaxt="n" , xlab="body mass (kg)" , ylab="brain volume (cc)" , col=rangi2 , pch=16 )
axis_unscale( 1 , at=quantile(d$mass) , d$mass )
axis_unscale( 2 , at=quantile(d$brain) , factor=max(d$brain) )

mass_seq <- seq(from=-1,to=1.5,length.out=30)
mu <- link(m7.2,data=list(mass_std=mass_seq))
mu <- apply(mu,2,mean)
lines( mass_seq , mu )

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

