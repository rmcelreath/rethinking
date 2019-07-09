\name{reedfrogs}
\alias{reedfrogs}
\docType{data}
\title{Data on reed frog predation experiments}
\description{
  Data on lab experiments on the density- and size-dependent
  predation rate of an African reed frog, \emph{Hyperolius
    spinigularis},
  from Vonesh and Bolker 2005 
}
\usage{data(reedfrogs)}
\format{
  Various data with variables:
  \describe{
    \item{\code{density}}{initial tadpole density (number of tadpoles
      in a 1.2 x 0.8 x 0.4 m tank) [experiment 1]}
    \item{\code{pred}}{factor: predators present or absent [experiment 1]}
    \item{\code{size}}{factor: big or small tadpoles [experiment 1]}
    \item{\code{surv}}{number surviving}
    \item{\code{propsurv}}{proportion surviving (=surv/density) [experiment 1]}
  }
}
\source{
  Vonesh and Bolker (2005) Compensatory larval responses shift
  trade-offs associated with predator-induced hatching plasticity.
  Ecology 86:1580-1591
}
\examples{
data(reedfrogs)
boxplot(propsurv~size*density*pred,data=reedfrogs)
}
\keyword{datasets}
