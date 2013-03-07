\name{coeftab}
\alias{coeftab}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Coefficient tables}
\description{
  Returns a table of model coefficients in rows and models in columns.
}
\usage{
coeftab( ... , se=FALSE , se.inside=FALSE , nobs=TRUE , digits=4 , rotate=FALSE , compare=FALSE )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{A series of fit models, separated by commas}
  \item{se}{Include standard errors in table?}
  \item{se.inside}{Print standard errors in same cell as estimates}
  \item{nobs}{Print number of observations for each model?}
  \item{digits}{Number of digits to round numbers to}
  \item{rotate}{If TRUE, rows are models and columns are coefficients}
  \item{compare}{Add AICc information?}
}
\details{
  This function provides a way to compare estimates across models.
}
\value{
}
\references{McElreath 2011, Statistical Rethinking.}
\author{Richard McElreath}
\seealso{}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

