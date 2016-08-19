\name{coerce_index}
\alias{coerce_index}
\alias{check_index}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Build and check integer index variables}
\description{
  These functions assist with building (\code{coerce_index}) and checking (\code{check_index}) integer index variables of the kind needed to define varying effect models.
}
\usage{
coerce_index( ... )
check_index( x )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{A comma-separted list of variables. See details.}
  \item{x}{A vector of integers to check for contiguity}
}
\details{
  Varying effect models often require index variables that begin at 1 and comprise only integers. These variables are used to lookup specific parameter values inside of the model. For example, it is common to define varying intercepts with a linear model like \code{a0 + a_id[id[i]]}. Here the variable \code{id} is an index variable. It has one value for each case, defining which individual applies to that case.

  When raw data exist as factors, these index variables much be converted to integers. This is trickier than it sounds, because R uses an internal integer represntation for factors, \code{levels}, that can conflict with ordinary integer representations.

  The function \code{coerce_index} deals with that complication. When the input is a single vector of factors, it returns an integer vector beginning at 1 and with contiguous values. 

  When the input is instead a comma-separated list of factors, it returns a list in which each factor has been converted to integers, but all levels in all factors were merged so that the same labels across factors always have the same integer values in the result. For example, suppose cases refer to dyads and there are two factors, \code{id1} and \code{id2}, that indicate which pair of individuals are present in each dyad. The labels in these variables should refer to the same individuals. Passing both simultaneously to \code{coerce_index} ensures that the results respect that fact.

  The function \code{check_index} merely checks an integer vector to see if it is contiguous.
}
\value{
For \code{coerce_index}, the result is either a single vector of integers (if the input was a single vector) or rather a list of vectors of integers (if the input was a list of vectors).
}
\references{}
\author{Richard McElreath}
\seealso{}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

