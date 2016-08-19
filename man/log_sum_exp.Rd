\name{log_sum_exp}
\alias{log_sum_exp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Rethinking utility functions}
\description{
  Numerically stable computation of log sums of exponentiated values.
}
\usage{
log_sum_exp( x )
}
\arguments{
    \item{x}{vector of values}
}
\details{
    This function is used by \code{WAIC} to compute the log average probability used in the formula for WAIC. In that context, a vector of log-probabilities is passed to \code{log_sum_exp}, obviating the need to explicitly exponentiate. This helps to avoid rounding errors that occur when working with direct probabilities.
}
%- maybe also 'usage' for other objects documented here.
\references{}
\author{Richard McElreath}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

