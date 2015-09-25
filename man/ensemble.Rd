\name{ensemble}
\alias{ensemble}
\title{Simulate an ensemble of posterior predictions}
\description{
  Uses \code{link} and \code{sim} for a list of \code{map} or \code{map2stan} model fits to construct Akaike weighted ensemble of predictions.
}
\usage{
ensemble( ... , data , n=1e3 , func=WAIC , weights , WAIC=TRUE , refresh=0 , 
          replace=list() , do_link=TRUE , do_sim=TRUE )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{\code{map} or \code{map2stan} models}
  \item{data}{Optional data to compute predictions over, as in \code{link} and \code{sim}}
  \item{n}{Number of samples to draw from posterior for each model}
  \item{func}{Function to use in computing criterion for model weights}
  \item{weights}{Optional vector of weights to use. If present, \code{func} is ignored.}
  \item{WAIC}{Deprecated: If \code{TRUE}, use \code{func} to compute weights. Otherwise tries to use DIC.}
  \item{refresh}{Progress update refresh interval. 0 suppresses output.}
  \item{replace}{Optional named list with replacement posterior samples. Used for maginalizing over random effects, for example. See example in \code{\link{link}}.}
  \item{do_link}{If \code{TRUE}, compute and return \code{link} results}
  \item{do_sim}{If \code{TRUE}, compute and return \code{sim} results}
}
\details{
  This function calls \code{\link{link}} and \code{\link{sim}} for each fit model given as input. The results are then combined into ensemble link and simulation output, where samples from each model are represented in proportion to the Akaike weights. Akaike weights are calculated by \code{\link{compare}}, using \code{func}, unless an explicit vector \code{weights} is provided. The values in \code{weights} will be normalized to sum to one, if they do not already.

  Note that no averaging is done by this function. Parameters are not averaged, and predictions are not averaged. 
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

