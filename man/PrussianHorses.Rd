\name{PrussianHorses}
\alias{PrussianHorses}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Numbers of soldiers killed by horse kicks in the Prussian army 1875-1894}
\description{
  Data collected by Ladsilaus von Bortkiewicz (1898), Das Gesetz der kleinen Zahlen. 
}
\usage{
data(PrussianHorses)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
}
\format{
    \enumerate{
        \item kicks: Count of soldier deaths from horse kicks
        \item year: Calendar year
        \item corps: Military corps (division). G is the guard corps.
    }
}
\value{
}
\references{von Bortkiewicz (1898), Das Gesetz der kleinen Zahlen}
\author{Richard McElreath}
\seealso{}
\examples{
data(PrussianHorses)
d <- PrussianHorses

dat <- list(
    kicks = d$kicks,
    year = as.integer(as.factor(d$year)), # year as category
    corps = as.integer(d$corps)
)

m0 <- ulam(
    alist(
        kicks ~ poisson( lambda ),
        log(lambda) <- d[corps] + y[year],
        d[corps] ~ normal(0,0.5),
        y[year] ~ normal(0,0.5)
    ), data=dat )

plot(precis(m0,2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

