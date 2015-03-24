\name{dordlogit}
\alias{dordlogit}\alias{rordlogit}\alias{pordlogit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Ordered categorical log-odds probability density}
\description{
  Functions for computing density, cumulative probability and producing random samples from an ordered categorical probability density.
}
\usage{
dordlogit( x , phi , a , log=FALSE )
pordlogit( x , phi , a , log=FALSE )
rordlogit( n , phi , a )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Integer values to compute densities or probabilies of}
  \item{a}{Vector of log-odds intercepts}
  \item{phi}{Linear model of log-odds}
  \item{log}{If \code{TRUE}, returns log-probability instead of probability}
}
\details{
  These functions provide for common density calculations for the ordered categorical probability density commonly known as an "ordered logit" or "ordered logistic." This is the same probability density assumed by the \code{polr} function (library \code{MASS}).
}
\references{}
\author{Richard McElreath}
\seealso{}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

