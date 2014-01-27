\name{map2stan}
\alias{map2stan}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Build RStan models from formulas}
\description{
  Compiles lists of formulas, like those used in \code{map}, into Stan model code. Allows for arbitary fixed effect and mixed effect regressions. Also computes DIC.
}
\usage{
map2stan( flist , data , start , sample=TRUE , iter=2000 , chains=1 , 
    debug=FALSE , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{flist}{A formula or list of formulas that define the likelihood and priors. See details.}
  \item{start}{A named list specifying (fixed) parameters and their initial values}
  \item{data}{A data frame or list containing the data}
  \item{sample}{If \code{FALSE}, builds Stan code without sampling}
  \item{iter}{Number of iterations of sampling. By default, half of these iterations are warmup}
  \item{chains}{Number of independent chains to sample from}
  \item{debug}{If \code{TRUE}, prints various internal steps to help with debugging}
  \item{...}{Additional arguments to pass to \code{\link{stan}}}
}
\details{
  This command provides a convenient interface for building arbitary fixed effect and mixed effect generalized linear models, as defined by a list of formulas. Syntax is similar to \code{map}, but also allowing multivariate priors corresponding to varying (aka random) effects.
  
  \code{flist} should be a either a single formula that defines the likelihood function or rather a list of formulas that define the likelihood and linear models and priors for parameters. See examples below.
  
  Likelihood formulas take the form \code{y ~ dfoo(bar)}, where \code{y} is the outcome variable, \code{dfoo} is a density function such as \code{dnorm}, and \code{bar} is a parameter of the density.
  
  Prior formulas take the same form, but the outcome should be a parameter name. Identical priors can be defined for multiple parameters by using \code{c(par1,par2,...)} on the left hand side of the formula. See example below.
  
  A special case of prior formula is for varying effects. For single varying effects, such as varying intercepts alone, all that is needed is to define a prior and mark it as conditional on a grouping variable in the data. For example: \code{aj|id ~ dnorm(0,sigma_id)} specifies a vector of varying effects \code{aj}, one for each unique value in \code{id}. For correlated varying effects, such as both varying intercepts and slopes, a parameter vector is specified and a multivariate prior is used instead. For example: \code{c(aj,bj)|id ~ dmvnorm(0,Sigma_id)} specifices varying intercepts \code{aj} and varying slopes \code{bj}
  
  Linear models can be specified as formulas of the form \code{mu ~ a + b*x} for a direct link. To use a link function, use the form \code{link(mu) ~ a + b*x} or \code{mu ~ invlink(a + b*x)}. The names "link" and "invlink" must be recognized by \code{map}. It currently recognizes \code{log}/\code{exp} and \code{logit}/\code{logistic}. Any other link function can be coded directly into the likelihood formula. For example \code{dfoo(par=invlink(mu))}.
  
  The Stan model code includes a generated quantities block that computes the deviance for each vector (iteration) of parameter samples. When sampling completes, \code{map2stan} computes DIC, the deviance information criterion, from the samples. DIC information is available from \code{show} and \code{DIC}, as well as being attributes of the returned object.
  
  Functions meant for processing \code{map2stan} fits are: \code{plotchains}, \code{extract.samples}, and \code{stancode}.
  
  Methods are defined for \code{coef}, \code{summary}, \code{logLik}, \code{vcov}, \code{nobs}, \code{deviance}, and \code{show}.
}
\value{
    Returns an object of class \code{map} with the following slots.
    \item{call}{The function call}
    \item{stanfit}{\code{stanfit} object returned by \code{\link{stan}}}
    \item{coef}{The posterior means}
    \item{vcov}{Minimal variance-covariance matrix, just holding diagonal variances}
    \item{data}{The data}
    \item{start}{List of starting values that was used in sampling}
    \item{formula}{Formula list from call}
    \item{formula_parsed}{List of parsed formula information. Useful mainly for debugging.}
}
\references{McElreath 2011, Statistical Rethinking.}
\author{Richard McElreath}
\seealso{\code{\link{map}},\code{\link{stan}}}
\examples{
library(rethinking)
data(chimpanzees)

# first a MAP fit of fixed effect model
# omitted prior for 'a' indicates improper flat prior
m1 <- map(
    list(
        pulled.left ~ dbinom(1,theta),
        logit(theta) ~ a + bp*prosoc.left + bpc*condition*prosoc.left,
        bp ~ dnorm(0,10),
        bpc ~ dnorm(0,10)
    ) ,
    data=chimpanzees,
    start=list(a=0,bp=0,bpc=0)
)

precis(m1)

# note that Stan doesn't allow "." in variable names
# we replace them with _
d <- list( 
    pulled_left = chimpanzees$pulled.left ,
    prosoc_left = chimpanzees$prosoc.left ,
    condition = chimpanzees$condition ,
    actor = as.integer( chimpanzees$actor )
)

# now RStan fit of same fixed effect model
m2 <- map2stan(
    list(
        pulled_left ~ dbinom(1,theta),
        logit(theta) ~ a + bp*prosoc_left + bpc*condition*prosoc_left ,
        bp ~ dnorm(0,10),
        bpc ~ dnorm(0,10)
    ) ,
    data=d,
    start=list(a=0,bp=0,bpc=0)
)

precis(m2)
plotchains(m2)
post2 <- extract.samples(m2)

# now RStan fit of model with varying intercepts on actor
m3 <- map2stan(
    list(
        pulled_left ~ dbinom(1,theta),
        logit(theta) ~ a + aj + bp*prosoc_left + bpc*condition*prosoc_left,
        aj|actor ~ dnorm( 0 , sigma_actor ),
        bp ~ dnorm(0,10),
        bpc ~ dnorm(0,10),
        sigma_actor ~ dcauchy(0,1)
    ) ,
    data=d,
    start=list(a=0,bp=0,bpc=0),
    iter=7000 , warmup=1000
)

precis(m3)
plotchains(m3)
post3 <- extract.samples(m3)

coeftab(m1,m2,m3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

