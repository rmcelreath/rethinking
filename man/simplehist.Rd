\name{simplehist}
\alias{simplehist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simple histograms}
\description{
  Simple integer-valued histograms, for displaying count distributions.
}
\usage{
simplehist( x , round=TRUE , ylab="Frequency" , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector of values to construct histogram from}
  \item{round}{When \code{TRUE}, rounds values in \code{x} before plotting}
  \item{ylab}{Label on vertical axis}
  \item{...}{Other parameters to pass to plot}
}
\details{
  This function constructs clean histograms for count data. Non-integer data can be rounded to nearest integer before plotting. Internally, this function is little more than \code{plot(table(x))}.
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

