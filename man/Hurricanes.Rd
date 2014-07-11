\name{Hurricanes}
\alias{Hurricanes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Hurricane fatalities and gender of names}
\description{
  Data used in Jung et al 2014 analysis of effect of gender of name on hurricane fatalities. Note that hurricanes Katrina (2005) and Audrey (1957) were removed from the data.
}
\usage{
data(Hurricanes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
}
\format{
    \enumerate{
        \item name : Given name of hurricane
        \item year : Year of hurricane
        \item deaths : number of deaths
        \item category : Severity code for storm
        \item min_pressure : Minimum pressure, a measure of storm strength; low is stronger
        \item damage_norm : Normalized estimate of damage in dollars
        \item female : Indicator variable for female name
        \item femininity : 1-11 scale from totally masculine (1) to totally feminine (11) for name. Average of 9 scores from 9 raters.
    }
}
\value{
}
\references{Jung et al. 2014. Female hurricanes are deadlier than male hurricanes. PNAS.}
\author{}
\seealso{}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

