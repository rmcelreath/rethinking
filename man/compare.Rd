\name{compare}
\alias{compare}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compare fit models using AICc/BIC}
\description{
  Returns a table of model comparison statistics
}
\usage{
compare( ... , nobs , sort="AICc" , delta=TRUE )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{A series of fit models, separated by commas}
  \item{nobs}{Number of observations to use in calculating AICc/BIC}
  \item{sort}{Which criterion to sort table by. "AICc" and "BIC" are valid.}
  \item{delta}{Include differences between AICc/BIC values}
}
\details{
  This function computes AICc and BIC values for fit models and returns a table sorted by ascending AICc. Each row in this table is a model, and the various columns provide AICc, BIC, numbers of parameters, and model weights.
}
\value{
}
\references{McElreath 2011, Statistical Rethinking.}
\author{Richard McElreath}
\seealso{}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

