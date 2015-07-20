\name{simplehist}
\alias{simplehist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simple histograms}
\description{
  Simple integer-valued histograms, for displaying count distributions.
}
\usage{
simplehist( x , ylab="Frequency" , xlab="Count" , ycounts=TRUE , adjust=1 , 
lcol="black" , bins=NULL , show.counts=0 , xlim=NULL , ylim=NULL , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector of values to construct histogram from}
  \item{ylab}{Label on vertical axis}
  \item{xlab}{Label on horizontal axis}
  \item{ycounts}{?}
  \item{lcol}{Line color}
  \item{bins}{Bins, defaults to zero and positive integers}
  \item{show.counts}{?}
  \item{...}{Other parameters to pass to plotting functions.}
}
\details{
  This function constructs clean histograms for count data. It needs to be rewritten, though.
}
\value{
}
\references{}
\author{Richard McElreath}
\seealso{\code{\link{hist}}}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

