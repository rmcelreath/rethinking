\name{dmvnorm2}
\alias{dmvnorm2}\alias{rmvnorm2}
\title{Multivariate Gaussian probability density}
\description{
  This is an alternative parameterization of the ordinary multivariate Gaussian probability density.
}
\usage{
dmvnorm2( x , Mu , sigma , Rho , log=FALSE )
rmvnorm2( n , Mu=rep(0,length(sigma)) , sigma=rep(1,length(Mu)) , 
          Rho=diag(length(Mu)) , method="chol" )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Values to compute densities of}
  \item{Mu}{Mean vector}
  \item{sigma}{Vector of standard deviations}
  \item{Rho}{Correlation matrix}
  \item{log}{If \code{TRUE}, returns log-density instead of density}
  \item{n}{Number of random observations to sample}
}
\details{
  These functions merely compose the variance-covariance matrix from separate standard deviation and correlation matrix arguments. They then use \code{dmvnorm} and \code{rmvnorm} from the \code{mvtnorm} package to perform calculations.
}
\references{}
\author{Richard McElreath}
\seealso{}
\examples{
dmvnorm2( c(1,0) , Mu=c(0,0) , sigma=c(1,1) , Rho=diag(2) )
rmvnorm2( 10 , Mu=c(1,2) )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

