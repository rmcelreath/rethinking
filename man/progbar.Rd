\name{progbar}
\alias{progbar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Progress display}
\description{
  Provides a progress display with estimated time remaining, assuming rate constant process.
}
\usage{
progbar(current, min = 0, max = 100, starttime,
                 update.interval = 100, show.rate = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{current}{Current loop index of process}
  \item{min}{Minimum loop index, usually 0 or 1}
  \item{max}{Maximum loop index}
  \item{starttime}{Time stamp when process began, from \code{Sys.time}}
  \item{update.interval}{How often to display progress}
}
\details{
  This function provides useful progress information and estimated time until completion for long looped operations, like sampling from MCMC. 
}
\value{
}
\references{}
\author{Richard McElreath}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

