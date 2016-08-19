\name{chimpanzees}
\alias{chimpanzees}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Chimpanzee prosocialty experiment data}
\description{
  Data from behavior trials in a captive group of chimpanzees, housed in Lousiana. From Silk et al. 2005. Nature 437:1357-1359.
}
\usage{
data(chimpanzees)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
}
\format{
    \enumerate{
        \item actor : name of actor
        \item recipient : name of recipient (NA for partner absent condition)
        \item condition : partner absent (0), partner present (1)
        \item block : block of trials (each actor x each recipient 1 time)
        \item trial : trial number (by chimp = ordinal sequence of trials for each chimp, ranges from 1-72; partner present trials were interspersed with partner absent trials)
        \item prosocial_left : 1 if prosocial (1/1) option was on left
        \item chose_prosoc : choice chimp made (0 = 1/0 option, 1 = 1/1 option)
        \item pulled_left : which side did chimp pull (1 = left, 0 = right)
    }
}
\value{
}
\references{Silk et al. 2005. Nature 437:1357-1359.}
\author{Richard McElreath}
\seealso{}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

