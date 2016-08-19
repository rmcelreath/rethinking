\name{shade}
\alias{shade}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Shade density intervals}
\description{
  Adds a shaded interval to an existing density plot or regression plot.
}
\usage{
shade(object, lim, label = NULL, col = col.alpha("black",
                 0.15), border = NA, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A \code{density} or \code{formula} object that defines the density OR a \code{matrix} object that defines the plot region. See details.}
  \item{lim}{For a density, a vector of two values, indicating upper and lower bounds of interval, on the horizontal axis. For a plot region, a list of x-axis values corresponding to y-axis points in the matrix object.}
  \item{label}{Optional label to center in interval.}
  \item{col}{Color to shade the interval. Default is transparent gray.}
  \item{border}{Border of shaded region. Default is no border, \code{NA}.}
  \item{...}{Other parameters to pass to \code{polygon}, which actually draws the shaded region.}
}
\details{
  This function uses \code{\link{polygon}} to draw a shaded region under a density curve or on a regression plot. The function assumes the plot already exists. See the examples below.
  
  There are two plotting interfaces, for densities. First, if the density is plotted from kernal estimation, using perhaps \code{\link{density}}, then the same density estimate should be passed as the first parameter. See example. Second, if the density is plotted from a series of x and y values, from perhaps a grid approximation, then a formula can be passed to define the curve. See example.
  
  For plotting confidence regions on a regression plot, the matrix object should contain y-axis values defining the border of the region. The first row of the matrix should be the bottom of the region, and the second row should be the top. Each column should correspond to the x-axis values in \code{lim}. See example.
}
\value{
}
\references{}
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

# plotting a shaded confidence interval on a regression plot
data(cars)
m <- lm( dist ~ speed , cars )
p <- extract.samples( m )
x.seq <- seq( from=min(cars$speed)-1 , to=max(cars$speed)+1 , length.out=30 )
mu.ci <- sapply( x.seq , function(x) PI( p[,1] + p[,2]*x ) )
plot( dist ~ speed , cars )
abline( m )
shade( mu.ci , x.seq )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

