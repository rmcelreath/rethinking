\name{mcreplicate}
\alias{mcreplicate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multi-core version of replicate}
\description{
  Uses the \code{parallel} library to distribute \code{\link{replicate}} processing across cores.
}
\usage{
mcreplicate(n, expr, refresh = 0.1, mc.cores=2 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Number of replications}
  \item{expr}{Code to replicate}
  \item{refresh}{Status update refresh interval}
  \item{mc.cores}{Number of cores to use}
}
\details{
  This function uses \code{\link{mclapply}} to distribute replications across cores. It then simplifies the result to an array.
}
\value{
}
\references{}
\author{Richard McElreath}
\seealso{\code{\link{mclapply}}, \code{\link{replicate}}}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

