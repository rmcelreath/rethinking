\name{KosterLeckie}
\alias{KosterLeckie}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Gift exchange data}
\description{
  Data on household gift exchanges from Koster and Leckie. There are two data frames: \code{kl_dyads} and \code{kl_households}.
}
\usage{
data(KosterLeckie)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
}
\format{
    \code{kl_dyads}:
    \enumerate{
        \item hidA : household ID for A in dyad
        \item hidB : household ID for B in dyad
        \item did : dyad ID number
        \item giftsAB : count of gifts from A to B
        \item giftsBA : count of gifts from B to A
        \item offset : relative rate of contact in dyad
        \item drel1 : 
        \item drel2 :
        \item drel3 :
        \item drel4 :
        \item dlndist :
        \item dass :
        \item d0125 :
    }
    \code{kl_households}:
    \enumerate{
        \item hid : household ID
        \item hgame : 
        \item hfish :
        \item hpigs :
        \item hwealth :
        \item hpastor :
    }
}
\value{
}
\references{Koster and Leckie (2014) Food sharing networks in lowland Nicaragua: An application of the social relations model to count data. Social Networks.}
\seealso{}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

