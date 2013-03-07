\name{shade}
\alias{shade}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Shade density intervals}
\description{
  Adds a shaded interval to an existing density plot.
}
\usage{
shade( object , lim , label=NULL , col="#00000066" , border=NA , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A \code{density} or \code{formula} object that defines the density. See details.}
  \item{lim}{A vector of two values, indicating upper and lower bounds of interval, on the horizontal axis.}
  \item{label}{Optional label to center in interval.}
  \item{col}{Color to shade the interval. Default is transparent gray.}
  \item{border}{Border of shaded region. Default is no border, \code{NA}.}
  \item{...}{Other parameters to pass to \code{polygon}, which actually draws the shaded region.}
}
\details{
  This function uses \code{\link{polygon}} to draw a shaded region under a density curve. The function assumes the density plot already exists. See the examples below.
  
  There are two plotting interfaces. First, if the density is plotted from kernal estimation, using perhaps \code{\link{density}}, then the same density estimate should be passed as the first parameter. See example. Second, if the density is plotted from a series of x and y values, from perhaps a grid approximation, then a formula can be passed to define the curve. See example.
}
\value{
}
\references{McElreath 2011, Statistical Rethinking.}
\author{Richard McElreath}
\seealso{\code{\link{density}}, \code{\link{dens}}}
\examples{
models <- seq( from=0 , to=1 , length.out=100 )
prior <- rep( 1 , 100 )
likelihood <- dbinom( 6 , size=9 , prob=models )
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)

# using formula interface
plot( posterior ~ models , type="l" )
shade( posterior ~ models , c(0,0.5) )

# using density interface
samples <- sample( models , size=1e4 , replace=TRUE , prob=posterior )
plot( density(samples) )
shade( density(samples) , PCI(samples) )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

