\name{drawdag}
\alias{drawdag}
\alias{drawopenpaths}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot Directed Acyclic Graphs (dagitty)}
\description{
  A fancier version of dagitty's plot function, as well as a way to show open paths given a conditioning set.
}
\usage{
drawdag( x , col_arrow="black" , col_segment="black" , col_labels="black" , cex=1 , 
  lwd=1.5 , goodarrow=TRUE , xlim , ylim , shapes , col_shapes , radius=3 , 
  add=FALSE , ... )

drawopenpaths( x , Z=list() , col_arrow="red" , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A dagitty graph}
  \item{col_arrow}{Color or vector of colors for the graph arrows}
  \item{col_segment}{Color or vector of colors for the graph segments}
  \item{col_labels}{Color or vector of colors for the graph text labels}
  \item{cex}{Size of text labels}
  \item{lwd}{Width of arrow lines}
  \item{goodarrow}{Use \code{Arrow} from \code{shape} package to draw arrows}
  \item{xlim}{Optional plot limits}
  \item{ylim}{Optional plot limits}
  \item{shapes}{A named list of variables with one of "c" for an open circle or "fc" for a filled circle}
  \item{col_shapes}{A named list of colors to correspond to the shapes list}
  \item{radius}{Radius of shapes circles}
  \item{add}{If TRUE, draw over existing DAG in active plot}
  \item{Z}{List of variables to condition on when computing open paths}
  \item{...}{Optional arguments to pass to other functions}
}
\details{
  \code{drawdag} is a modified version of \code{plot.dagitty} but with additional stylistic options. By default, it draws arrows in black and with thicker line width. It also uses the nicer arrows drawn by \code{Arrows} in the \code{shape} package. If the DAG doesn't already have \code{\link{coordinates}}, then \code{\link{graphLayout}} is called to provide them.

  \code{drawopenpaths} uses \code{\link{dagitty::paths}} to compute and then overdraw open paths, given an exposure and outcome variable. It uses \code{drawdag} to perform the drawing. Requires that the DAG already be displayed. Only open paths are drawn in the overlay color.
}
\value{
}
\references{}
\author{Richard McElreath}
\seealso{}
\examples{
\dontrun{

ex1 <- dagitty("dag {
    X [exposure]
    Y [outcome]
    U [unobserved]
    Z -> X -> Y
    X <- U -> Y
}")
coordinates(ex1) <- list( x=c(Z=0,X=1,Y=1,U=0) , y=c(Z=0,U=0.5,X=0,Y=1) )
drawdag( ex1 )

# example of drawing open paths
drawdag( ex1 )
drawopenpaths( ex1 ) # open backdoor

drawdag( ex1 , col_arrow="gray" )
drawopenpaths( ex1 , Z=list("U") ) # closed backdoor
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

