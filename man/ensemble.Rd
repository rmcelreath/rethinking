\name{ensemble}
\alias{ensemble}
\title{Simulate an ensemble of posterior predictions}
\description{
  Uses \code{link} and \code{sim} for a list of \code{map} or \code{map2stan} model fits to construct Akaike weighted ensemble of predictions.
}
\usage{
ensemble( ... , data , n=1e3 , WAIC=TRUE , refresh=0 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{\code{map} or \code{map2stan} models}
  \item{data}{Optional data to compute predictions over, as in \code{link} and \code{sim}}
  \item{n}{Number of samples to draw from posterior for each model}
  \item{WAIC}{Use WAIC to compute weights. Otherwise tries to use DIC.}
  \item{refresh}{Progress update refresh interval. 0 suppresses output.}
}
\details{
  This function calls \code{\link{link}} and \code{\link{sim}} for each fit model given as input. The results are then combined into ensemble link and simulation output, where samples from each model are represented according the Akaike weights.
}
\value{
}
\references{
}
\author{Richard McElreath}
\seealso{\code{\link{link}}, \code{\link{sim}}, \code{\link{compare}}}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

