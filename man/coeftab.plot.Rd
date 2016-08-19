\name{coeftab_plot}
\alias{coeftab_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plots of coefficient tables}
\description{
  Plots coefficient tables produced by \code{coeftab}, clustered either by models or by parameter names.
}
\usage{
coeftab_plot( x , y , pars , col.ci="black" , by.model=FALSE , 
    prob=0.95 , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Object produced by \code{coeftab}}
  \item{y}{\code{NULL} and unused. Required for compatibility with base \code{plot}}
  \item{pars}{Optional vector of parameter names or indexes to display. If missing, all parameters shown.}
  \item{col.ci}{Color to draw confidence intervals}
  \item{by.model}{Cluster estimates by model instead of by parameter (default)}
  \item{prob}{Probability mass for confidence intervals. Default is 0.95.}
}
\details{
  This function plots the tabular output of \code{\link{coeftab}}, using a \code{\link{dotchart}}. By default, estimates are grouped by parameter, with a row for each model. Model's without a parameter still appear as a row, but with no estimate. By setting \code{by.model=TRUE}, the dotchart will instead be grouped by model, with each row being a parameter. 

  MAP estimates are displayed with percentile confidence (credible) intervals. Default is 95\% intervals. Use \code{prob} to change the interval mass.
}
\value{
}
\references{}
\author{Richard McElreath}
\seealso{\code{\link{coeftab}}, \code{\link{dotchart}}}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

