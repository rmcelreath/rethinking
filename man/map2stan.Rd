\name{map2stan}
\alias{map2stan}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Build RStan models from formulas}
\description{
  Compiles lists of formulas, like those used in \code{map}, into Stan model code. Allows for arbitary fixed effect and mixed effect regressions. Computes DIC and WAIC. Allows for simple imputation of missing values.
}
\usage{
map2stan( flist , data , start , pars , constraints=list() , types=list() , 
  sample=TRUE , iter=2000 , warmup=floor(iter/2) , chains=1 , debug=FALSE , 
  verbose=FALSE , WAIC=TRUE , cores=1 , rng_seed , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{flist}{A formula or list of formulas that define the likelihood and priors. Can also pass in a \code{map} model fit. See details.}
  \item{data}{A data frame or list containing the data}
  \item{start}{Optional named list specifying parameters and their initial values}
  \item{pars}{Optional: character vector of parameters to return samples for}
  \item{constraints}{Optional: named list of custom parameter constraints, using Stan notation}
  \item{types}{Optional: named list of custom parameter types, using Stan notation}
  \item{sample}{If \code{FALSE}, builds Stan code without sampling}
  \item{iter}{Number of iterations of sampling. By default, half of these iterations are warmup.}
  \item{warmup}{Number of warmup iterations. By default, half of \code{iter}.}
  \item{chains}{Number of independent chains to sample from}
  \item{debug}{If \code{TRUE}, prints various internal steps to help with debugging}
  \item{verbose}{If \code{TRUE}, prints extra progress messages.}
  \item{WAIC}{When \code{TRUE}, computes WAIC after sampling, storing the result}
  \item{cores}{Number of processor cores to distribute chains over, using \code{\link{parallel}}.}
  \item{...}{Additional arguments to pass to \code{\link{stan}}}
}
\details{
  This command provides a convenient interface for building arbitary fixed effect and mixed effect generalized linear models, as defined by a list of formulas. Syntax is similar to \code{map}, but also allowing multivariate priors corresponding to varying (aka random) effects, as well as simple imputation schemes.
  
  \code{flist} should be either (1) a single formula that defines the likelihood function or rather a list of formulas that define the likelihood and linear models and priors for parameters (see examples below) or (2) a previously fit \code{map} model.
  
  Likelihood formulas take the form \code{y ~ dfoo(bar)}, where \code{y} is the outcome variable, \code{dfoo} is a density function such as \code{dnorm}, and \code{bar} is a parameter of the density.
  
  Prior formulas take the same form, but the outcome should be a parameter name. Identical priors can be defined for multiple parameters by using \code{c(par1,par2,...)} on the left hand side of the formula. See example below.
  
  A special case of prior formula is for varying effects. For single varying effects, such as varying intercepts alone, all that is needed is to define a prior and mark it as conditional on a grouping variable in the data. For example: \code{aj[id] ~ dnorm(0,sigma_id)} specifies a vector of varying effects \code{aj}, one for each unique value in \code{id}. For correlated varying effects, such as both varying intercepts and slopes, a parameter vector is specified and a multivariate prior is used instead. For example: \code{c(aj,bj)[id] ~ dmvnorm(0,Sigma_id)} specifices varying intercepts \code{aj} and varying slopes \code{bj}.
  
  Linear models can be specified as formulas of the form \code{mu <- a + b*x} for a direct link. To use a link function, use the form \code{link(mu) <- a + b*x}. The name "link" must be recognized by \code{map2stan}. It currently recognizes \code{log} and \code{logit}.
  
  Imputation of missing values is available by specifying distributions for predictor variables that contain \code{NA} values. \code{map2stan} will split the variable into observed values and a vector of parameters used to estimate the missing values, assigning the same distribution to each. See the example.
  
  The \code{start} list is optional. When missing from the list, for each parameter with a defined prior, an initial value will be sampled from the prior. Sampled initial values will try to respect parameter constraints. For varying effect parameter vectors, initial values will always be set to zero. Specific initial values can be specified in the \code{start} list. See examples below. 
  
  The Stan model code includes a generated quantities block that computes the deviance for each iteration of parameter samples. When sampling completes, \code{map2stan} computes DIC, the deviance information criterion, from the samples. DIC information is available from \code{show} and \code{DIC}, as well as being attributes of the returned object. 
  
  WAIC can be computed with \code{WAIC}, or by setting \code{WAIC=TRUE} when calling \code{map2stan}. This is currently the default. WAIC is calculated entirely after Stan completes sampling. 
  
  Methods are defined for \code{\link{extract.samples}}, \code{\link{link}}, \code{\link{sim}}, \code{\link{ensemble}}, \code{\link{compare}}, \code{coef}, \code{summary}, \code{logLik}, \code{vcov}, \code{nobs}, \code{deviance}, \code{plot}, \code{pairs}, and \code{show}.
}
\value{
    Returns an object of class \code{map2stan} with the following slots.
    \item{call}{The function call}
    \item{model}{Stan model code}
    \item{stanfit}{\code{stanfit} object returned by \code{\link{stan}}}
    \item{coef}{The posterior means}
    \item{vcov}{Minimal variance-covariance matrix, just holding diagonal variances}
    \item{data}{The data}
    \item{start}{List of starting values that were used in sampling}
    \item{pars}{Parameter names monitored in samples}
    \item{formula}{Formula list from call}
    \item{formula_parsed}{List of parsed formula information. Useful mainly for debugging.}
}
\references{}
\author{Richard McElreath}
\seealso{\code{\link{resample}}, \code{\link{map}}, \code{\link{stan}}, \code{\link{link}}, \code{\link{sim}}, \code{\link{glimmer}}}
\examples{
\dontrun{
library(rethinking)
data(chimpanzees)

# don't want any variables with NAs
d <- list( 
    pulled_left = chimpanzees$pulled_left ,
    prosoc_left = chimpanzees$prosoc_left ,
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
    data=d, chains=2, cores=1 )

precis(m2)
summary(m2)
plot(m2)
pairs(m2)

# now RStan fit of model with varying intercepts on actor
m3 <- map2stan(
    alist(
        pulled_left ~ dbinom(1,theta),
        logit(theta) <- a + aj[actor] + bp*prosoc_left + bpc*condition*prosoc_left,
        aj[actor] ~ dnorm( 0 , sigma_actor ),
        a ~ dnorm(0,10),
        bp ~ dnorm(0,10),
        bpc ~ dnorm(0,10),
        sigma_actor ~ dcauchy(0,1)
    ) ,
    data=d,
    iter=5000 , warmup=1000 , chains=2 , cores=1 )

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
    iter=5000 , warmup=1000 , chains=2 , cores=1 )

precis(m4)
summary(m4)
plot(m4)

# compare posterior means
coeftab(m2,m3,m4)
plot(coeftab(m2,m3,m4))

# show WAIC for m2,m3,m4
compare(m2,m3,m4)
plot(compare(m2,m3,m4))

###########
# varying slopes models

# varying slopes on actor
# also demonstrates use of multiple linear models
# see Chapter 13 for discussion
m5 <- map2stan(
    alist(
        # likeliood
        pulled_left ~ dbinom(1,p),

        # linear models
        logit(p) <- A + (BP + BPC*condition)*prosoc_left,
        A <- a + a_actor[actor],
        BP <- bp + bp_actor[actor],
        BPC <- bpc + bpc_actor[actor],

        # adaptive prior
        c(a_actor,bp_actor,bpc_actor)[actor] ~
                                dmvnorm2(0,sigma_actor,Rho_actor),

        # fixed priors
        c(a,bp,bpc) ~ dnorm(0,1),
        sigma_actor ~ dcauchy(0,2),
        Rho_actor ~ dlkjcorr(4)
    ) , data=d , iter=5000 , warmup=1000 , chains=3 , cores=3 )

# same model but with non-centered parameterization
# see Chapter 13 for explanation and more elaborate example

m6 <- map2stan(
    alist(
        # likeliood
        pulled_left ~ dbinom(1,p),

        # linear models
        logit(p) <- A + (BP + BPC*condition)*prosoc_left,
        A <- a + a_actor[actor],
        BP <- bp + bp_actor[actor],
        BPC <- bpc + bpc_actor[actor],

        # adaptive prior - non-centered
        c(a_actor,bp_actor,bpc_actor)[actor] ~
                                dmvnormNC(sigma_actor,Rho_actor),

        # fixed priors
        c(a,bp,bpc) ~ dnorm(0,1),
        sigma_actor ~ dcauchy(0,2),
        Rho_actor ~ dlkjcorr(4)
    ) , data=d , iter=5000 , warmup=1000 , chains=3 , cores=3 )

###########
# Imputation example

# simulate data:
#  linear regression with two predictors
#  both predictors have valules missing at random
N <- 100
N_miss <- 10
x1 <- rnorm( N )
x2 <- rnorm( N )
y <- rnorm( N , 2*x1 - 0.5*x2 , 1 )
x1[ sample(1:N,size=N_miss) ] <- NA
x2[ sample(1:N,size=N_miss) ] <- NA

# formula with distributions assigned to both predictors
f <- alist(
    y ~ dnorm( mu , sigma ),
    mu <- a + b1*x1 + b2*x2,
    x1 ~ dnorm( mu_x1, sigma_x1 ),
    x2 ~ dnorm( mu_x2, sigma_x2 ),
    a ~ dnorm( 0 , 100 ),
    c(b1,b2) ~ dnorm( 0  , 10 ),
    c(mu_x1,mu_x2) ~ dnorm( 0 , 100 ),
    c(sigma_x1,sigma_x2) ~ dcauchy(0,2),
    sigma ~ dcauchy(0,2)
)

m <- map2stan( f , data=list(y=y,x1=x1,x2=x2) , sample=TRUE )

# show observed outcomes against retrodicted outcomes
# cases with missing values shown with red posterior intervals
v <- link(m)
mu <- apply( v , 2 , mean )
ci <- apply( v , 2 , PI )
plot( y ~ mu )
cicols <- ifelse( is.na(x1) | is.na(x2) , "red" , "gray" )
for( i in 1:N ) lines( ci[,i] , rep(y[i],2) , col=cicols[i] )
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

