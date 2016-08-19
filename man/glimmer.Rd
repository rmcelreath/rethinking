\name{glimmer}
\alias{glimmer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{glm/glmer formulas to map/map2stan formulas}
\description{
  Converts a \code{glm} or \code{glmer} model formula into a formula list of the kind used by \code{map} and \code{map2stan}
}
\usage{
glimmer( formula , data , family=gaussian , prefix=c("b_","v_") , 
    default_prior="dnorm(0,10)" , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{A formula of the sort used by \code{\link{glm}} or \code{glmer}}
  \item{data}{A data frame or list containing the data}
  \item{family}{Name of family, as in \code{glm}/\code{glmer}. Recognizes \code{gaussian}, \code{binomial}, and \code{poisson} with arbitrary links}
  \item{prefix}{Text prefixes for fixed and varying effects parameters}
  \item{default_prior}{Text of default prior for regression parameters}
  \item{...}{Additional parameters to pass to \code{map2stan}}
}
\details{
This function parses and converts a \code{glmer}-style mixed model formula into a list of formulas appropriate for \code{map2stan} input. Interactions are pre-multiplied, and factors are automatically converted into a series of dummy variables. In the absence of varying effects in the input formula, the resulting formula list will also be suitable for \code{map} input in most cases.

The function \code{glmer} is provided by package \code{lme4}, but it is not required for \code{glimmer} to work.
}
\value{
    Returns an invisible list with two slots:
    \item{f}{list of formulas, suitable for passing to \code{map2stan}}
    \item{d}{list of data, suitable for passing to \code{map2stan}}
}
\references{}
\author{Richard McElreath}
\seealso{\code{\link{map2stan}}}
\examples{
\dontrun{
library(rethinking)
data(UCBadmit)

# varying intercepts
f3 <- cbind(admit,reject) ~ (1|dept) + applicant.gender
m3 <- glimmer( f3 , UCBadmit , binomial )
m3s <- map2stan( m3$f , data=m3$d )

# varying intercepts and slopes
f4 <- cbind(admit,reject) ~ (1+applicant.gender|dept) + applicant.gender
m4 <- glimmer( f4 , UCBadmit , binomial )
m4s <- map2stan( m4$f , data=m4$d )
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

