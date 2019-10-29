\name{cv_quap}
\alias{cv_quap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Diagnostic trace and rank histogram plots for MCMC output}
\description{
  Cross validation for \code{quap} model fits.
}
\usage{
cv_quap( quap_model, lno = 1, pw = FALSE, cores = 1, ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{quap_model}{\code{quap} model fit}
  \item{lno}{Number of observations to leave out in each fold}
  \item{pw}{Pointwise results (TRUE) or summed across folds (FALSE)}
  \item{cores}{Number of cores to use. If great than 1, uses \code{mclappy} with folds}
  \item{...}{Additional arguments to pass to \code{\link{quap}}}
}
\details{
  This function constructs cross validation folds from an existing \code{quap} model fit and associated data. It then fits the model to each fold and returns either the fit to each fold (when \code{pw=TRUE}) or the summed performance across all folds.

  The default is leave-one-out cross-validation, when \code{lno=1}.
}
\references{}
\author{Richard McElreath}
\seealso{}
\examples{
\dontrun{

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

