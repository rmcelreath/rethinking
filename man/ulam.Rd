\name{ulam}
\alias{ulam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Build RStan models from formulas}
\description{
  Compiles lists of formulas into Stan model code. Allows for arbitary fixed effect and mixed effect regressions. Allows for explicit typing of variables within the formula list. Much more flexible than \code{map2stan}.
}
\usage{
ulam( flist , data , pars , pars_omit , start , chains=1 , cores=1 , iter=1000 , 
  control=list(adapt_delta=0.95) , distribution_library=ulam_dists , 
  macro_library=ulam_macros , custom , declare_all_data=TRUE , log_lik=FALSE , 
  sample=TRUE , messages=TRUE , pre_scan_data=TRUE , coerce_int=TRUE , 
  sample_prior=FALSE , file=NULL , cmdstan=FALSE , threads=1 ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{flist}{A formula or list of formulas that define the likelihood and priors. Can also pass in a \code{quap} or a previous \code{ulam} model fit. See details.}
  \item{data}{A data frame or list containing the data}
  \item{pars}{Optional: character vector of parameters to return samples for}
  \item{pars_omit}{Optional: character vector of parameters to exclude from samples}
  \item{start}{Optional. Either: (1) a named list specifying parameters and their initial values or (2) a function to return such a named list.}
  \item{chains}{Number of independent chains to sample from}
  \item{cores}{Number of processor cores to distribute chains over.}
  \item{iter}{Number of iterations of sampling. By default, half of these iterations are warmup.}
  \item{control}{Optional list of control parameters for \code{stan}. Default increases target acceptance rate (\code{adapt_delta}) to 0.95.}
  \item{distribution_library}{List of distribution templates.}
  \item{macro_library}{List of function and distribution macros.}
  \item{custom}{Optional list of custom resources. See details and examples.}
  \item{declare_all_data}{When \code{TRUE}, all variables in the data list are declared in the Stan model code. When \code{FALSE}, only used variables are declared.}
  \item{log_lik}{Return log likelihood of each observation in samples. Used for calculating WAIC and LOO.}
  \item{sample}{If \code{FALSE}, builds Stan code without sampling}
  \item{messages}{Show various warnings and informational messages}
  \item{pre_scan_data}{Scan data at start to (1) check for variables that are integer but not type integer and (2) strip any scale() attributes}
  \item{coerce_int}{If \code{pre_scan_data}, forces integer variables to be type integer}
  \item{sample_prior}{If \code{TRUE}, removes data probabilities from model to sample only prior distribution of parameters. Used by corresponding \code{extract.prior} method.}
  \item{messages}{If \code{TRUE}, prints various internal steps to help with debugging}
  \item{file}{If character string X, loads X.rds instead of fitting when exists; otherwise saves result to X.rds}
  \item{cmdstan}{When TRUE, uses cmdstanr instead of rstan to run chains. To make cmdstan the default engine, use \code{set_ulam_cmdstan(TRUE)}.}
  \item{threads}{When threads > 1, attempts to multithread individual chains using Stan's reduce_sum function. Requires cmdstan=TRUE.}
  \item{...}{Additional arguments to pass to \code{\link{stan}}}
}
\details{
  \code{ulam} provides a slim version of Stan syntax that omits blocks but still allows for explicit variable types and dimensions. The basic model formula is the same as \code{map2stan}, but the syntax does not assume a GLMM structure. This allows it to be more flexible, but it also means more mistakes are possible. With great power comes great responsibility. 

  The function of a model formula is to related the variables to one another. There are three types of variables: (1) data, (2) parameters, and (3) local variables. Model formulas are composed of multiple lines. Each line defines a variable using either a distributional assumption like:

  \code{y ~ normal( mu , sigma )}

  or rather a determinsitic assignment like:

  \code{mu <- a + b*x}

  The basic structure of such a definition is:

  \code{type[dimension]: name[dimension] ~ distribution( arguments )}

  The \code{type} delcaration is optional. So the most basic defintion can be just \code{y ~ bernoulli(theta)}, but a full declariation can be more detailed, when necessary. The examples below show how matrix variables can be defined in this syntax.

  For determinstic assignments like:

  \code{mu <- a + b*x}

  It is also possible to use a control word to specify how the values are returned. Using \code{save>} returns the values in the posterior samples. For example:

  \code{save> mu <- a + b*x}

  will return \code{mu} for each case and posterior sample. This works by duplicating the code in both the model block, where it is used to compute the log-probability, and in generated quantities.

  It is also possible to use \code{gq} to evaluate the assignment only after sampling, in Stan's generated quantities block. This is useful for derived values that are not needed in computing the posterior but may be useful afterwards. For example, constrasts could be calculated this way. In the examples, the line:

  \code{gq> bp_diff <- bp[1] - bp[2]}

  is used to calculate the posterior distribution of the difference between the two parameters. The code is added to Stan's generated quantities, so that it doesn't slow down the model block.

  The control tag \code{transpars} can be used to place an assignment in Stan's transformed parameters block. Keep in mind that any other intermediate calculations must also be placed in the same block. Finally, \code{transdata} places the assignment in the transformed data block. This means it will only execute once, before sampling begins. But it also means that the values will not be available post-sampling for helper functions like \code{link}. As such, it will usually be better to transform data before passing it into \code{ulam}.

  When \code{cmdstan=TRUE}, the \code{cmdstanr} package will be used instead \code{rstan} to compile and sample from models. This is generally superior, as more recent versions of Stan can be used this way. But some features are not yet implemented, such as passing custom inits to the chains. To make cmdstan the default engine, use \code{set_ulam_cmdstan(TRUE)}. You can then ignore the \code{cmdstan} argument when calling \code{ulam}.

  The use of \code{cmdstan=TRUE} is also the only way to currently use multi-threading of individual chains, using the \code{threads} argument. When \code{cmdstan=TRUE} and \code{threads} is set greater than 1, \code{ulam} will try to recode the model so that each chain is spread over multiple cores. This can easily halve sampling time. At the moment, this only works for models with a single outcome variable. 
  
  Methods are defined for \code{\link{extract.samples}}, \code{\link{extract.prior}}, \code{\link{link}}, \code{\link{sim}}, \code{\link{compare}}, \code{coef}, \code{summary}, \code{logLik}, \code{lppd}, \code{vcov}, \code{nobs}, \code{deviance}, \code{WAIC}, \code{PSIS}, \code{plot}, \code{traceplot}, \code{trankplot}, \code{pairs}, and \code{show}.
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
# also recode condition to an index {1,0} -> {1,2}
d <- list( 
    pulled_left = chimpanzees$pulled_left ,
    prosoc_left = chimpanzees$prosoc_left ,
    condition = as.integer( 2 - chimpanzees$condition ) ,
    actor = as.integer( chimpanzees$actor ) ,
    blockid = as.integer( chimpanzees$block )
)

# simple logistic regression
m1 <- ulam(
    alist(
        pulled_left ~ bernoulli(theta),
        logit(theta) <- a + bp[condition]*prosoc_left  ,
        a ~ normal(0,4),
        bp[condition] ~ normal(0,1)
    ) ,
    data=d, chains=2, cores=1 , sample=TRUE )

precis(m1,depth=2)
plot(m1,depth=2)
pairs(m1)

# same model, but save theta so it is return in samples
# note 'save>' in second line of formula
m1b <- ulam(
    alist(
        pulled_left ~ bernoulli(theta),
        save> logit(theta) <- a + bp[condition]*prosoc_left  ,
        a ~ normal(0,4),
        bp[condition] ~ normal(0,1)
    ) ,
    data=d, chains=2, cores=1 , sample=TRUE )

# same model, but use gq to compute contrast between conditions
# note that order does matter. bp_diff should come before bp[] is defined
m1c <- ulam(
    alist(
        pulled_left ~ bernoulli(theta),
        logit(theta) <- a + bp[condition]*prosoc_left  ,
        gq> bp_diff <- bp[1] - bp[2],
        a ~ normal(0,4),
        bp[condition] ~ normal(0,1)
    ) ,
    data=d, chains=2, cores=1 , sample=TRUE )

# can also transform data inside model, using transdata> tag.
# this is more efficient, because it only evaluates once, not during sampling.
# for example, this constructs prosoc_right variable:
m1d <- ulam(
    alist(
        pulled_left ~ bernoulli(theta),
        logit(theta) <- a + bp[condition]*prosoc_right  ,
        transdata> prosoc_right <- 1 - prosoc_left,
        a ~ normal(0,4),
        bp[condition] ~ normal(0,1)
    ) ,
    data=d, chains=2, cores=1 , sample=TRUE )

# now model with varying intercepts on actor
m2 <- ulam(
    alist(
        pulled_left ~ bernoulli(theta),
        logit(theta) <- a + aj[actor] + bp[condition]*prosoc_left,
        aj[actor] ~ normal( 0 , sigma_actor ),
        a ~ normal(0,4),
        bp[condition] ~ normal(0,1),
        sigma_actor ~ exponential(1)
    ) ,
    data=d, chains=2 , cores=1 , sample=TRUE )

precis(m2)
plot(m2)

# varying intercepts on actor and experimental block
m3 <- ulam(
    alist(
        pulled_left ~ bernoulli(theta),
        logit(theta) <- a + aj[actor] + ak[blockid] + bp[condition]*prosoc_left,
        aj[actor] ~ normal( 0 , sigma_actor ),
        ak[blockid] ~ normal( 0 , sigma_block ),
        a ~ dnorm(0,4),
        bp[condition] ~ dnorm(0,1),
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
        pulled_left ~ bernoulli(theta),

        # linear models
        logit(theta) <- A + BP*prosoc_left,
        A <- a + v[actor,1],
        BP <- bp + v[actor,condition+1],

        # adaptive prior
        vector[3]: v[actor] ~ multi_normal( 0 , Rho_actor , sigma_actor ),

        # fixed priors
        c(a,bp) ~ normal(0,1),
        sigma_actor ~ exponential(1),
        Rho_actor ~ lkjcorr(4)
    ) , data=d , chains=3 , cores=1 , sample=TRUE )

# same model but with non-centered parameterization
# see Chapter 13 for explanation and more elaborate example
m4 <- ulam(
    alist(
        # likeliood
        pulled_left ~ bernoulli(theta),

        # linear models
        logit(theta) <- A + BP*prosoc_left,
        A <- a + v[actor,1],
        BP <- bp + v[actor,condition+1],

        # adaptive prior
        matrix[actor,3]: v <- compose_noncentered( sigma_actor , L_Rho_actor , z ),
        matrix[3,actor]: z ~ normal( 0 , 1 ),

        # fixed priors
        c(a,bp) ~ normal(0,1),
        vector[3]: sigma_actor ~ exponential(1),
        cholesky_factor_corr[3]: L_Rho_actor ~ lkj_corr_cholesky( 4 )
    ) , data=d , chains=3 , cores=1 , sample=TRUE )

# same as m5, but without hiding the construction of v
m5 <- ulam(
    alist(
        # likeliood
        pulled_left ~ bernoulli(theta),

        # linear models
        logit(theta) <- A + BP*prosoc_left,
        A <- a + v[actor,1],
        BP <- bp + v[actor,condition+1],

        # adaptive prior
        matrix[actor,3]: v <- t(diag_pre_multiply( sigma_actor , L_Rho_actor ) * z),
        matrix[3,actor]: z ~ normal( 0 , 1 ),

        # fixed priors
        c(a,bp,bpc) ~ normal(0,1),
        vector[3]: sigma_actor ~ exponential(1),
        cholesky_factor_corr[3]: L_Rho_actor ~ lkj_corr_cholesky( 4 )
    ) , data=d , chains=3 , cores=1 , sample=TRUE )

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

