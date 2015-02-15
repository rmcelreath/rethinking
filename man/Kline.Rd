\name{Kline}
\alias{Kline}\alias{Kline2}\alias{islandsDistMatrix}\alias{islands}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Oceanic islands data}
\description{
  Data on historic tool complexity and demography in various Oceanic islands societies. There are three data tables: (1) \code{Kline} is the basic data table; (2) \code{Kline2} contains latitude and longitude columns; and (3) \code{islandsDistMatrix} is a matrix of pairwise distances between islands, in thousands of kilometers.
}
\usage{
data(Kline)
data(Kline2)
data(islandsDistMatrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
}
\format{
    \enumerate{
        \item culture : Name of island culture
        \item population : Historical population size
        \item contact : low or high contact rate with other islands
        \item total_tools : number of tools in historical tool kit
        \item mean_TU : another measure of tool complexity
        \item lat : latitude of island
        \item lon : longitude of island
        \item lon2 : longitude coded for ease of plotting
        \item logpop : natural logarithm of population
    }
}
\value{
}
\references{Kline, M.A. and R. Boyd. 2010. Proc R Soc B 277:2559--2564.}
\seealso{}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

