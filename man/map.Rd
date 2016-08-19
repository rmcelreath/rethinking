\name{map}
\alias{map}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find maximum a posteriori (MAP) estimates}
\description{
  Find mode of posterior distribution for arbitrary fixed effect models.
}
\usage{
map( flist , data , start , method="BFGS" , hessian=TRUE , debug=FALSE , verbose=FALSE , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{flist}{A formula or \code{alist} of formulas that define the likelihood and priors. See details.}
  \item{data}{A data frame or list containing the data}
  \item{start}{Optional named list specifying parameters and their initial values}
  \item{method}{Search method for \code{optim}. Defaults to \code{BFGS}.}
  \item{hessian}{If \code{TRUE}, attempts to compute the Hessian}
  \item{debug}{If \code{TRUE}, prints various internal steps to help with debugging}
  \item{...}{Additional arguments to pass to \code{optim}}
}
\details{
  This command provides a convenient interface for finding maximum a posteriori estimates for models defined by a formula or by a list of formulas. This procedure is equivalent to penalized maximum likelihood estimation, and therefore can be used to define many common regularization procedures.
  
  \code{flist} should be a either a single formula that defines the likelihood function or rather a list of formulas that define the likelihood and priors for parameters. See examples below.
  
  Likelihood formulas take the form \code{y ~ dfoo(bar)}, where \code{y} is the outcome variable, \code{dfoo} is a density function such as \code{dnorm}, and \code{bar} is a parameter of the density.
  
  Prior formulas take the same form, but the outcome should be a parameter name. Identical priors can be defined for multiple parameters by using \code{c(par1,par2,...)} on the left hand side of the formula. See example below.
  
  Linear models can be specified as formulas of the form \code{mu <- a + b*x} for a direct link. To use a link function, use the form \code{link(mu) <- a + b*x} or \code{mu <- invlink(a + b*x)}. The names "link" and "invlink" must be recognized by \code{map}. It currently recognizes \code{log}/\code{exp} and \code{logit}/\code{logistic}. Any other link function can be coded directly into the likelihood formula. For example \code{y ~ dfoo(invlink(mu))}.
  
  The \code{start} list is optional. For each parameter with a defined prior, an initial value will be sampled from the prior. Explicit named initial values can also be provided in this list.
  
  Methods are defined for \code{coef}, \code{summary}, \code{logLik}, \code{vcov}, \code{nobs}, \code{deviance}, \code{link}, \link{DIC}, and \code{show}.
}
\value{
    Returns an object of class \code{map} with the following slots.
    \item{call}{The function call}
    \item{coef}{The estimated posterior modes}
    \item{vcov}{Variance-covariance matrix}
    \item{optim}{List returned by \code{optim}}
    \item{data}{The data}
    \item{formula}{Formula list}
    \item{fminuslogl}{Minus log-likelihood function that accepts a vector of parameter values}
}
\references{}
\author{Richard McElreath}
\seealso{\code{\link{optim}}}
\examples{
data(cars)
flist <- alist(
    dist ~ dnorm( mean=mu , sd=sigma ) ,
    mu <- a+b*speed ,
    c(a,b) ~ dnorm(0,10) , 
    sigma ~ dcauchy(0,1)
)
fit <- map( flist , start=list(a=40,b=0.1,sigma=20) , data=cars )

# regularized logistic regression example
y <- c( rep(0,10) , rep(1,10) )
x <- c( rep(-1,9) , rep(1,11) )

flist0 <- alist(
    y ~ dbinom( prob=p , size=1 ) ,
    logit(p) <- a + b*x
)

flist1 <- alist(
    y ~ dbinom( prob=p , size=1 ),
    logit(p) <- a + b*x ,
    c(a,b) ~ dnorm(0,10)
)

# without priors, same problem as:
# glm3a <- glm( y ~ x , family=binomial , data=list(y=y,x=x) )
fit3a <- map( flist0 , data=list(y=y,x=x) , start=list(a=0,b=0) )
precis(fit3a)

# now with regularization
fit3b <- map( flist1 , data=list(y=y,x=x) , start=list(a=0,b=0) )
precis(fit3b)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

