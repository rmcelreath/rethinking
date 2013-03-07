\name{glmm}
\alias{glmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generalized linear multilevel models}
\description{
  Convenient interface to fitting generalized linear multilevel (mixed) models with package \code{lme4}.
}
\usage{
glmm( ... , REML=FALSE )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{Arguments to pass on to \code{glmer}}
  \item{REML}{If \code{TRUE}, fits a Gaussian glmm using restricted maximum likelihood}
}
\details{
  This function merely provides a convenient interface for fitting models with \code{glmer} from the \code{lme4} package. All \code{glmm} really does is make \code{REML=FALSE} the default.
}
\value{
}
\references{McElreath 2011, Statistical Rethinking.}
\author{Richard McElreath}
\seealso{\code{\link{lmer}}}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

