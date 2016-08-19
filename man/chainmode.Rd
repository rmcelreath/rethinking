\name{chainmode}
\alias{chainmode}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find mode of a continuous density estimate}
\description{
  Returns estimated mode of a density computed from samples.
}
\usage{
chainmode( chain , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chain}{Values, e.g. sampled from a posterior via MCMC}
  \item{...}{Optional arguments passed to density calculation}
}
\details{
  This function just finds the x value that maximizes the y density in the density estimate.
}
\value{
}
\references{}
\author{Richard McElreath}
\seealso{}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

