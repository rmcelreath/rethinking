\name{plot.map2stan}
\alias{plot.map2stan}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Trace plot of map2stan chains}
\description{
  Shows trace plots for Stan samples produced by a \code{map2stan} fit, annotated by number of effective samples. Automatic paging and adjustable number of columns and rows per page.
}
\usage{
plot( object , pars , col=rethink_palette , alpha=1 , bg=gray(0.7,0.5) , ask=TRUE , window , n_cols=3 , max_rows=5 , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{map2stan model fit}
  \item{pars}{Character vector of parameters to display}
  \item{col}{Vector of colors to use for chains. Recycled.}
  \item{alpha}{alpha transparency of chains}
  \item{bg}{Background fill for warmup samples}
  \item{ask}{Whether to pause for paging. Default is \code{TRUE}.}
  \item{window}{Range of samples to display. Default is all samples.}
  \item{n_cols}{Number of columns per page}
  \item{max_rows}{Maximum number of rows per page}
  \item{...}{Additional arguments to pass to plot commands}
}
\details{
  This is the default trace plot method for \code{\link{map2stan}} model fits. 
}
\value{
}
\references{}
\author{Richard McElreath}
\seealso{}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

