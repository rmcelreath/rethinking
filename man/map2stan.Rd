\name{map2stan}
\alias{map2stan}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Build RStan models from formulas}
\description{
  Compiles lists of formulas, like those used in \code{map}, into Stan model code. Allows for arbitary fixed effect and mixed effect regressions. Also computes DIC and WAIC.
}
\usage{
map2stan( flist , data , start , pars , constraints=list() , types=list() , 
    sample=TRUE , iter=2000 , chains=1 , debug=FALSE , WAIC=FALSE , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{flist}{A formula or list of formulas that define the likelihood and priors. See details.}
  \item{data}{A data frame or list containing the data}
  \item{start}{A named list specifying parameters and their initial values}
  \item{pars}{Optional: character vector of parameters to return samples for}
  \item{constraints}{Optional: named list of custom parameter constraints, using Stan notation}
  \item{types}{Optional: named list of custom parameter types, using Stan notation}
  \item{sample}{If \code{FALSE}, builds Stan code without sampling}
  \item{iter}{Number of iterations of sampling. By default, half of these iterations are warmup}
  \item{chains}{Number of independent chains to sample from}
  \item{debug}{If \code{TRUE}, prints various internal steps to help with debugging}
  \item{WAIC}{When \code{TRUE}, computes WAIC after sampling, storing the result}
  \item{...}{Additional arguments to pass to \code{\link{stan}}}
}
\details{
  This command provides a convenient interface for building arbitary fixed effect and mixed effect generalized linear models, as defined by a list of formulas. Syntax is similar to \code{map}, but also allowing multivariate priors corresponding to varying (aka random) effects.
  
  \code{flist} should be a either a single formula that defines the likelihood function or rather a list of formulas that define the likelihood and linear models and priors for parameters. See examples below.
  
  Likelihood formulas take the form \code{y ~ dfoo(bar)}, where \code{y} is the outcome variable, \code{dfoo} is a density function such as \code{dnorm}, and \code{bar} is a parameter of the density.
  
  Prior formulas take the same form, but the outcome should be a parameter name. Identical priors can be defined for multiple parameters by using \code{c(par1,par2,...)} on the left hand side of the formula. See example below.
  
  A special case of prior formula is for varying effects. For single varying effects, such as varying intercepts alone, all that is needed is to define a prior and mark it as conditional on a grouping variable in the data. For example: \code{aj|id ~ dnorm(0,sigma_id)} specifies a vector of varying effects \code{aj}, one for each unique value in \code{id}. For correlated varying effects, such as both varying intercepts and slopes, a parameter vector is specified and a multivariate prior is used instead. For example: \code{c(aj,bj)|id ~ dmvnorm(0,Sigma_id)} specifices varying intercepts \code{aj} and varying slopes \code{bj}.
  
  Linear models can be specified as formulas of the form \code{mu ~ a + b*x} for a direct link. To use a link function, use the form \code{link(mu) ~ a + b*x}. The name "link" must be recognized by \code{map2stan}. It currently recognizes \code{log} and \code{logit}.
  
  The Stan model code includes a generated quantities block that computes the deviance for each vector (iteration) of parameter samples. When sampling completes, \code{map2stan} computes DIC, the deviance information criterion, from the samples. DIC information is available from \code{show} and \code{DIC}, as well as being attributes of the returned object.
  
  Functions meant for (re)processing \code{map2stan} fits are: \code{extract.samples}, \code{\link{resample}}, and \code{stancode}.
  
  Methods are defined for \code{coef}, \code{summary}, \code{logLik}, \code{vcov}, \code{nobs}, \code{deviance}, \code{plot}, \code{pairs}, and \code{show}.
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
\seealso{\code{\link{resample}},\code{\link{map}},\code{\link{stan}}}
\examples{
library(rethinking)
data(chimpanzees)

# note that Stan doesn't allow "." in variable names
# we replace them with _
# also don't want any variables with NAs
d <- list( 
    pulled_left = chimpanzees$pulled.left ,
    prosoc_left = chimpanzees$prosoc.left ,
    condition = chimpanzees$condition ,
    actor = as.integer( chimpanzees$actor ) ,
    blockid = as.integer( chimpanzees$block )
)

# RStan fit
m2 <- map2stan(
    alist(
        pulled_left ~ dbinom(1,theta),
        logit(theta) <- a + bp*prosoc_left + bpc*condition*prosoc_left ,
        a ~ dnorm(0,10),
        bp ~ dnorm(0,10),
        bpc ~ dnorm(0,10)
    ) ,
    data=d,
    start=list(a=0,bp=0,bpc=0)
)

precis(m2)
summary(m2)
plot(m2)
pairs(m2)

# now RStan fit of model with varying intercepts on actor
# note initial values for each varying intercept in start
m3 <- map2stan(
    alist(
        pulled_left ~ dbinom(1,theta),
        logit(theta) <- a + aj + bp*prosoc_left + bpc*condition*prosoc_left,
        aj[actor] ~ dnorm( 0 , sigma_actor ),
        a ~ dnorm(0,10),
        bp ~ dnorm(0,10),
        bpc ~ dnorm(0,10),
        sigma_actor ~ dcauchy(0,1)
    ) ,
    data=d,
    start=list(a=0,bp=0,bpc=0,sigma_actor=1,aj=rep(0,max(d$actor))),
    iter=7000 , warmup=1000 , chains=2
)

precis(m3)
plot(m3)
pairs(m3)

# varying intercepts on actor and experimental block
m4 <- map2stan(
    alist(
        pulled_left ~ dbinom(1,theta),
        logit(theta) <- a + aj + ak + bp*prosoc_left + bpc*condition*prosoc_left,
        aj[actor] ~ dnorm( 0 , sigma_actor ),
        ak[blockid] ~ dnorm( 0 , sigma_block ),
        a ~ dnorm(0,10),
        bp ~ dnorm(0,10),
        bpc ~ dnorm(0,10),
        sigma_actor ~ dcauchy(0,1),
        sigma_block ~ dcauchy(0,1)
    ) ,
    data=d,
    start=list(a=0,bp=0,bpc=0,sigma_actor=1,sigma_block=1,aj=rep(0,7),ak=rep(0,max(d$blockid))),
    iter=20000 , warmup=5000 , chains=2
)

precis(m4)
summary(m4)
plot(m4)

# compare posterior means
coeftab(m2,m3,m4)

# show DIC for m2,m3,m4
sapply( list(m2,m3,m4) , DIC )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

