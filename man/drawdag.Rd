\name{drawdag}
\alias{drawdag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot Directed Acyclic Graph (dagitty)}
\description{
  A fancier version of dagitty's plot function.
}
\usage{
drawdag( x , col_arrow="black" , col_segment="black" , col_labels="black"  , 
    cex=1 , lwd=1.5 , goodarrow=TRUE , xlim , ylim , ... )
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
  \item{...}{Optional arguments to pass to other functions}
}
\details{
  This is a modified version of \code{plot.dagitty} but with additional stylistic options. By default, it draws arrows in black and with thicker line width. It also uses the nicer arrows drawn by \code{Arrows} in the \code{shape} package.
}
\value{
}
\references{}
\author{Richard McElreath}
\seealso{}
\examples{
\dontrun{
library(dagitty)
library(shape)
ex1 <- dagitty("dag {
    Z -> X -> Y
    X <- U -> Y
}")
coordinates(ex1) <- list( x=c(Z=0,X=1,Y=1,U=0) , y=c(Z=0,U=0.5,X=0,Y=1) )
drawdag( ex1 , col_arrow=c("black","black","black","red") )
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

