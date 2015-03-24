\name{contour_xyz}
\alias{contour_xyz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Contour plot from equal length x,y,z vectors}
\description{
  Provides an interface to use \code{contour} by providing three equal length vectors for x, y and z coordinates.
}
\usage{
contour_xyz( x , y , z , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector of x values}
  \item{y}{vector of y values}
  \item{z}{vector of z values}
  \item{...}{other parameters to pass to \code{contour}}
}
\details{
  This function merely constructs a matrix suitable for \code{contour}, using x, y and z coordinates.
}
\value{
}
\references{}
\author{Richard McElreath}
\seealso{\code{\link{contour}}}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

