\name{link}
\alias{link}
\title{Intermediate predictions for map and map2stan models}
\description{
  Computes inverse-link linear model values for \code{map} and \code{map2stan} samples.
}
\usage{
link( fit , data , n=1000 , probs=NULL , refresh=0.1 , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{Object of class \code{map2stan}}
  \item{data}{Optional list of data to compute predictions over}
  \item{n}{Number of samples to use}
  \item{probs}{Not yet implemented}
  \item{refresh}{Refresh interval for progress display. Set to \code{refresh=0} to suppress display.}
  \item{...}{Other parameters to pass to someone}
}
\details{
  This function computes the value of each linear model at each sample for each case in the data. Inverse link functions are applied, so that for example a logit link linear model produces probabilities, using the logistic. 

  This function is used internally by \code{\link{WAIC}}. It's still a work in progress. Most notably, it doesn't work quite right for ordered logistic models.
  
  The \code{rethinking} package defines a generic function \code{link}, so methods can be defined for other model fit classes. I might eventually build methods for \code{lm} and \code{glm}.
}
\value{
}
\references{
}
\author{Richard McElreath}
\seealso{\code{\link{map}},\code{\link{map2stan}}}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

