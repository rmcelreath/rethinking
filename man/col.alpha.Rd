\name{col.alpha}
\alias{col.alpha}\alias{col.desat}\alias{col.dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Color utility functions}
\description{
  Functions for calculating transparent and desaturated colors.
}
\usage{
col.alpha( acol , alpha = 0.2 )
col.desat( acol , amt = 0.5 )
col.dist( x , mu = 0 , sd = 1 , col="slateblue" )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{acol}{A color name or RGB color}
  \item{alpha}{alpha transparency, where 1 means fully opaque and 0 fully transparent}
  \item{amt}{amount of desaturation of color to apply, where 1 means totally desaturated (grayscale)}
  \item{x}{A vector of values to use for calculating distances. See details below.}
  \item{mu}{Value (or vector of values) to use for calculating distances}
  \item{sd}{Standard deviation of distance function used to calculate transparency}
  \item{col}{A color to apply transparency to, based on distance}
}
\details{
  These functions allow for calculating transparency and desaturation for colors. \code{col.alpha} makes a base color transparent, while \code{col.desat} makes a color have less saturation.
  
  \code{col.dist} is used to make a list of transparent colors of differing alpha level. The levels are chosen based upon Gaussian distance from a chosen value, \code{mu}, with a chosen width of the function that determines how quickly colors become fully transparent, \code{sd}. For example, if \code{x} contains a column of data, then \code{col.dist} will return a vector of same length with transparency increasing as each value in \code{x} is distant from \code{mu}. This is useful for plotting data that emphasize points near some value or values.
}
\references{}
\author{Richard McElreath}
\seealso{\code{\link{col2rgb}}, \code{\link{rgb2hsv}}, \code{\link{rgb}}}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

