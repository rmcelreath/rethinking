\name{UFClefties}
\alias{UFClefties}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Ultimate Fighting Championship matches}
\description{
  Outcomes of televised Ultimate Fighting Championship (UFC) matches, as a function of handedness of fighters. Data coded by Pollet et al (2013).
}
\usage{
data(UFClefties)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
}
\format{
    \enumerate{
        \item fight : Unique identifier for match
        \item episode : Identifier for UFC episode
        \item fight.in.episiode : Order of fight in episode
        \item fighter1.win : 1 if fighter 1 won the match; 0 if fight 2 won
        \item fighter1 : Unique identifier for fighter 1
        \item fighter2 : Unique identifier for fighter 2
        \item fighter1.lefty : 1 if fighter 1 was left handed; 0 otherwise
        \item fighter2.lefty : 1 if fighter 2 was left handed; 0 otherwise
    }
}
\value{
}
\references{Pollet et al. 2013. Animal Behaviour 86:839--843.}
\seealso{}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

