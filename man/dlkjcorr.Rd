\name{dlkjcorr}
\alias{dlkjcorr}\alias{rlkjcorr}
\title{LKJ correlation matrix probability density}
\description{
  Functions for computing density and producing random samples from the LKJ onion method correlation matrix distribution.
}
\usage{
dlkjcorr( x , eta=1 , log=TRUE )
rlkjcorr( n , K , eta=1 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Matrix to compute probability density for}
  \item{eta}{Parameter controlling shape of distribution}
  \item{K}{Dimension of correlation matrix}
  \item{log}{If \code{TRUE}, returns log-probability instead of probability}
  \item{n}{Number of random matrices to sample}
}
\details{
  The LKJ correlation matrix distribution is based upon work by Lewandowski, Kurowicka, and Joe. When the parameter \code{eta} is equal to 1, it defines a flat distribution of correlation matrices. When \code{eta > 1}, the distribution is instead concentrated towards to identity matrix. When \code{eta < 1}, the distribution is more concentrated towards extreme correlations at -1 or +1.
  
  It can be easier to understand this distribution if we recognize that the individual correlations within the matrix follow a beta distribution defined on -1 to +1. Thus \code{eta} resembles \code{theta} in the beta parameterization with a mean p and scale (sample size) theta. 
  
  The function \code{rlkjcorr} returns an 3-dimensional array in which the first dimension indexes matrices. In the event that \code{n=1}, it returns instead a single matrix.
}
\references{Lewandowski, Kurowicka, and Joe. 2009. Generating random correlation matrices based on vines and extended onion method. Journal of Multivariate Analysis. 100:1989-2001.

Stan Modeling Language User's Guide and Reference Manual, v2.6.2}
\author{Richard McElreath}
\seealso{}
\examples{
R <- rlkjcorr(n=1,K=2,eta=4)
dlkjcorr(R,4)

# plot density of correlation
R <- rlkjcorr(1e4,K=2,eta=4)
dens( R[,1,2] )

# visualize 3x3 matrix
R <- rlkjcorr(1e3,K=3,eta=2)
plot( R[,1,2] , R[,1,3] , col=col.alpha("black",0.2) , pch=16 )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

