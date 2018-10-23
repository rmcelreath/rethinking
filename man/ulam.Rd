\name{ulam}
\alias{ulam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Build RStan models from formulas}
\description{
  Compiles lists of formulas into Stan model code. Allows for arbitary fixed effect and mixed effect regressions. Allows for explicit typing of variables within the formula list. Much more flexible than \code{map2stan}.
}
\usage{
ulam( flist , data , pars , start , chains=1 , cores=1 , iter=1000 , 
  control=list(adapt_delta=0.95) , distribution_library=ulam_dists , 
  macro_library=ulam_macros , declare_all_data=TRUE ,
  log_lik=FALSE , sample=TRUE , messages=TRUE , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{flist}{A formula or list of formulas that define the likelihood and priors. Can also pass in a \code{map} model fit. See details.}
  \item{data}{A data frame or list containing the data}
  \item{pars}{Optional: character vector of parameters to return samples for}
  \item{start}{Optional named list specifying parameters and their initial values}
  \item{chains}{Number of independent chains to sample from}
  \item{cores}{Number of processor cores to distribute chains over.}
  \item{iter}{Number of iterations of sampling. By default, half of these iterations are warmup.}
  \item{control}{Optional list of control parameters for \code{stan}. Default increases target acceptance rate (\code{adapt_delta}) to 0.95.}
  \item{distribution_library}{List of distribution templates.}
  \item{macro_library}{List of function and distribution macros.}
  \item{declare_all_data}{When \code{TRUE}, all variables in the data list are declared in the Stan model code. When \code{FALSE}, only used variables are declared.}
  \item{log_lik}{Return log likelihood of each observation in samples. Used for calculating WAIC and LOO.}
  \item{sample}{If \code{FALSE}, builds Stan code without sampling}
  \item{messages}{If \code{TRUE}, prints various internal steps to help with debugging}
  \item{WAIC}{When \code{TRUE}, computes WAIC after sampling, storing the result}
  \item{constraints}{Optional: named list of custom parameter constraints, using Stan notation}
  \item{...}{Additional arguments to pass to \code{\link{stan}}}
}
\details{
  \code{ulam} provides a slim version of Stan syntax that omits blocks but still allows for explicit variable types and dimensions. The basic model formula is the same as \code{map2stan}, but the syntax does not assume a GLMM structure. This allows it to be more flexible, but it also means more mistakes are possible. With great power comes great responsibility. 
  
  Methods are defined for \code{\link{extract.samples}}, \code{\link{extract.prior}}, \code{\link{link}}, \code{\link{sim}}, \code{\link{ensemble}}, \code{\link{compare}}, \code{coef}, \code{summary}, \code{logLik}, \code{vcov}, \code{nobs}, \code{deviance}, \code{plot}, \code{traceplot}, \code{pairs}, and \code{show}.
}
\value{
    Returns an object of class \code{ulam} with the following slots.
    \item{call}{The function call}
    \item{model}{Stan model code}
    \item{stanfit}{\code{stanfit} object returned by \code{\link{stan}}}
    \item{coef}{The posterior means}
    \item{vcov}{k-by-1 matrix containing the variance of each of k variables in posterior}
    \item{data}{The data}
    \item{start}{List of starting values that were used in sampling}
    \item{pars}{Parameter names monitored in samples}
    \item{formula}{Formula list from call}
    \item{formula_parsed}{List of parsed formula information. Useful mainly for debugging. Needed by helper functions.}
}
\references{}
\author{Richard McElreath}
\seealso{\code{\link{quap}}, \code{\link{map2stan}}, \code{\link{stan}}}
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

# simple logistic regression
m1 <- ulam(
    alist(
        pulled_left ~ binomial(1,theta),
        logit(theta) <- a + bp*prosoc_left + bpc*condition*prosoc_left ,
        a ~ normal(0,4),
        bp ~ normal(0,1),
        bpc ~ normal(0,1)
    ) ,
    data=d, chains=2, cores=1 , sample=TRUE )

precis(m1)
plot(m1)
pairs(m1)

# now model with varying intercepts on actor
m2 <- ulam(
    alist(
        pulled_left ~ binomial(1,theta),
        logit(theta) <- a + aj[actor] + bp*prosoc_left + bpc*condition*prosoc_left,
        aj[actor] ~ normal( 0 , sigma_actor ),
        a ~ normal(0,4),
        bp ~ normal(0,1),
        bpc ~ normal(0,1),
        sigma_actor ~ exponential(1)
    ) ,
    data=d, chains=2 , cores=1 , sample=TRUE )

precis(m2)
plot(m2)

# varying intercepts on actor and experimental block
m3 <- ulam(
    alist(
        pulled_left ~ binomial(1,theta),
        logit(theta) <- a + aj[actor] + ak[blockid] + bp*prosoc_left + bpc*condition*prosoc_left,
        aj[actor] ~ normal( 0 , sigma_actor ),
        ak[blockid] ~ normal( 0 , sigma_block ),
        a ~ dnorm(0,4),
        bp ~ dnorm(0,1),
        bpc ~ dnorm(0,1),
        sigma_actor ~ exponential(1),
        sigma_block ~ exponential(1)
    ) ,
    data=d, chains=2 , cores=1 , sample=TRUE )

precis(m3)
summary(m3)
plot(m3)

###########
# varying slopes models

# varying slopes on actor
# also demonstrates use of multiple linear models
# see Chapter 13 for discussion
m3 <- ulam(
    alist(
        # likeliood
        pulled_left ~ binomial(1,p),

        # linear models
        logit(p) <- A + (BP + BPC*condition)*prosoc_left,
        A <- a + a_actor[actor],
        BP <- bp + bp_actor[actor],
        BPC <- bpc + bpc_actor[actor],

        # adaptive prior
        c(a_actor,bp_actor,bpc_actor)[actor] ~
                                multi_normal(0,Rho_actor,sigma_actor),

        # fixed priors
        c(a,bp,bpc) ~ normal(0,1),
        sigma_actor ~ exponential(1),
        Rho_actor ~ lkjcorr(4)
    ) , data=d , chains=3 , cores=1 , sample=TRUE )

# same model but now with explicit vector of varying effects for each actor
# note indexing in linear models
m4 <- ulam(
    alist(
        # likeliood
        pulled_left ~ binomial(1,p),

        # linear models
        logit(p) <- A + (BP + BPC*condition)*prosoc_left,
        A <- a + v[actor,1],
        BP <- bp + v[actor,2],
        BPC <- bpc + v[actor,3],

        # adaptive prior
        vector[3]:v[actor] ~ multi_normal(0,Rho_actor,sigma_actor),

        # fixed priors
        c(a,bp,bpc) ~ normal(0,1),
        sigma_actor ~ exponential(1),
        Rho_actor ~ lkjcorr(4)
    ) , data=d , chains=3 , cores=1 , sample=TRUE )

# same model but with non-centered parameterization
# see Chapter 13 for explanation and more elaborate example
m5 <- ulam(
    alist(
        # likeliood
        pulled_left ~ binomial(1,p),

        # linear models
        logit(p) <- A + (BP + BPC*condition)*prosoc_left,
        A <- a + v[actor,1],
        BP <- bp + v[actor,2],
        BPC <- bpc + v[actor,3],

        # adaptive prior
        matrix[actor,3]: v <- compose_noncentered( sigma_actor , L_Rho_actor , z ),
        matrix[3,actor]: z ~ normal( 0 , 1 ),

        # fixed priors
        c(a,bp,bpc) ~ normal(0,1),
        vector[3]:sigma_actor ~ exponential(1),
        cholesky_factor_corr[3]: L_Rho_actor ~ lkj_corr_cholesky( 4 )
    ) , data=d , chains=3 , cores=1 , sample=TRUE )

# same as m5, but without hiding the construction of v
m5b <- ulam(
    alist(
        # likeliood
        pulled_left ~ binomial(1,p),

        # linear models
        logit(p) <- A + (BP + BPC*condition)*prosoc_left,
        A <- a + v[actor,1],
        BP <- bp + v[actor,2],
        BPC <- bpc + v[actor,3],

        # adaptive prior
        matrix[actor,3]: v <- t(diag_pre_multiply( sigma_actor , L_Rho_actor ) * z),
        matrix[3,actor]: z ~ normal( 0 , 1 ),

        # fixed priors
        c(a,bp,bpc) ~ normal(0,1),
        vector[3]:sigma_actor ~ exponential(1),
        cholesky_factor_corr[3]: L_Rho_actor ~ lkj_corr_cholesky( 4 )
    ) , data=d , chains=3 , cores=1 , sample=TRUE )

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

