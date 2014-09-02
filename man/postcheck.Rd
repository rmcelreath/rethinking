\name{postcheck}
\alias{postcheck}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot posterior predictive check}
\description{
  Plots a series of graphical posterior predictive checks for a \code{map} or \code{map2stan} model fit.
}
\usage{
postcheck( fit , prob=0.9 , window=20 , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{A \code{map} or \code{map2stan} model fit}
  \item{prob}{Width of prediction interval to display}
  \item{window}{Number of cases to display in each plot (for paging)}
  \item{...}{Additional arguments to pass someplace}
}
\details{
  This function uses \code{\link{link}} and \code{\link{sim}} internally to simulate posterior predictions for a fit model. These are then plotted over the observed outcomes used in fitting.
}
\value{
}
\references{McElreath 2011, Statistical Rethinking.}
\author{Richard McElreath}
\seealso{\code{\link{link}},\code{\link{sim}}}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

