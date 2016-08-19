\name{image_xyz}
\alias{image_xyz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Heat map from equal length x,y,z vectors}
\description{
  Provides an interface to use \code{image} by providing three equal length vectors for x, y and z coordinates.
}
\usage{
image_xyz( x , y , z , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector of x values}
  \item{y}{vector of y values}
  \item{z}{vector of z values}
  \item{...}{other parameters to pass to \code{image}}
}
\details{
  This function merely constructs a matrix suitable for \code{image}, using x, y and z coordinates.
}
\value{
}
\references{}
\author{Richard McElreath}
\seealso{\code{\link{image}}}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

