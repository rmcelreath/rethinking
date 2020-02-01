\name{Boxes}
\alias{Boxes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Social learning experimental data}
\description{
  Data from (and Stan models for) a cross-cultural experiment investigating the development of social learning in children
}
\usage{
data(Boxes)
data(Boxes_model)
data(Boxes_model_age)
data(Boxes_model_gender)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
}
\format{
    \enumerate{
        \item y : Outcome, one of 1=unchosen option, 2=majority option, or 3=minority option
        \item gender : Index of child's gender. 1=girl. 2=boy.
        \item age : Age of each child in years
        \item majority_first : whether the majority option was demonstrated first
        \item culture : ID of the site, from 1 to 8
    }
}
\value{
}
\references{van Leeuwen et al 2018. The development of human social learning across seven societies.  DOI: 10.1038/s41467-018-04468-2.}
\author{Richard McElreath}
\seealso{}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

