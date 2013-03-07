\name{precis}
\alias{precis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Precis of model fit}
\description{
  Displays concise parameter estimate information for an existing model fit.
}
\usage{
precis( model , type.s=FALSE , ci=TRUE , level=0.95 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{Fit model object}
  \item{type.s}{Show Type-S error probabilities for each estimate}
  \item{ci}{Show quadratic estimate confidence intervals}
  \item{level}{Width of confidence intervals}
}
\details{
  Creates a table of estimates and standard errors, with optional confidence intervals and Type-S error probabilities. Confidence intervals are quadratic estimates, derived from standard errors. Type-S probabilities are also derived from normal approximation to posterior of each parameter.
  
  Supported model classes include \code{lm}, \code{mle2}, \code{mer}, and \code{polr}. Can also provide expected value, standard deviation, and HPDI columns for a data frame, on the assumption that the columns in the data frame are the product of a Markov chain.
}
\value{
  A data frame with a row for each parameter.
}
\references{McElreath 2011, Statistical Rethinking.}
\author{Richard McElreath}
\seealso{\code{\link{summary}}}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

