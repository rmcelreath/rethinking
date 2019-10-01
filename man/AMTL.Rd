\name{AMTL}
\alias{AMTL}
\alias{AMTL_short}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Ante-mortem Tooth Loss Data}
\description{
  Data from four primate genera on tooth loss and its relationship to age and sex. Used for measurement error example in the textbook.
}
\usage{
data(AMTL_short)
data(AMTL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
}
\format{
    \enumerate{
        \item tooth_class : One of Anterior, Posterior, or Premolar
        \item specimen : Unique identifier for specimen
        \item num_amtl : Number of teeth missing of given class
        \item sockets : number of observable sockets that could be scored for missing teeth
        \item age : Estimated age of speciment at death
        \item stdev_age : Assigned uncertainty of age at death
        \item prob_male : Estimate of sex of specimen
        \item genus : Specimen genus, one of Homo, Pan, Papio, or Pongo
        \item population : Region specimen originates from
    }
}
\value{
}
\references{Gilmore, C.C. 2013. A Comparison of Antemortem Tooth Loss in Human Hunter-Gatherers and Non-human Catarrhines: Implications for the Identification of Behavioral Evolution in the Fossil Record. American Journal of Physical Anthropology. DOI: 10.1002/ajpa.22275.}
\seealso{}
\examples{
data(AMTL)
# plot proportion lost against age
plot( d$num_amtl / d$sockets , d$age )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

