\name{ulam}
\alias{ulam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Build RStan models from formulas}
\description{
  Compiles lists of formulas into Stan model code. Allows for arbitary fixed effect and mixed effect regressions. Much more flexible than \code{map2stan}.
}
\usage{
ulam( flist , data , pars , start , chains=1 , iter=1000 , 
  control=list(adapt_delta=0.95) , distribution_library=ulam_dists , 
  macro_library=ulam_macros , declare_all_data=TRUE ,
  log_lik=FALSE , sample=FALSE , messages=TRUE , ... )
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
  \item{rng_seed}{Optional explicit seed.}
  \item{rawstanfit}{When \code{TRUE}, returns just the \code{stanfit} object, as if the model had been called using \code{\link{stan}}.}
  \item{control}{Optional list of control parameters for \code{stan}. Default increases target acceptance rate (\code{adapt_delta}) to 0.95.}
  \item{add_unique_tag}{When \code{TRUE}, adds a comment to the Stan model code with the date-time stamp. This makes each model unique and will force Stan to recompile. Useful for avoiding segfault bugs when reusing compiled objects.}
  \item{code}{Optional list of custom Stan code to insert in model. See details and example.}
  \item{log_lik}{Return log likelihood of each observation in samples. Used for calculating WAIC and LOO.}
  \item{DIC}{Return deviance and DIC. This is deprecated and may be removed in a later version.}
  \item{declare_all_data}{When \code{TRUE}, all variables in the data list are declared in the Stan model code. When \code{FALSE}, only used variables are declared.}
  \item{do_discrete_imputation}{When \code{TRUE}, samples for missing binary predictors are returned. Not necessary to marginalize over discrete missing values.}
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

  When predictor variables are binary (0/1), \code{map2stan} will attempt to marginalize over any missing values. This is accomplished by building a mixture likelihood. Missingness in more than one binary variable can be accommodated this way, by building a list of all combinations of missingness among the variables and then a correspond vector of likelihood terms. The resulting Stan code contains a loop that computes the proper mixture and adds it to the target with \code{log_sum_exp}. The user may need to use the optional \code{constraints} list to bound hyperparameters. See the example.
  
  The \code{start} list is optional. When missing from the list, for each parameter with a defined prior, an initial value will be sampled from the prior. Sampled initial values will try to respect parameter constraints. For varying effect parameter vectors, initial values will always be set to zero. Specific initial values can be specified in the \code{start} list. See examples below. 

  The optional \code{code} argument can be used to pass a list of custom Stan code to be inserted into specific parts of the model. The list should be a list of lists. Each list should have the format \code{list("code",block="block",section="section",pos="pos")}. The first argument is the code to insert, as a character string. The named \code{block} slot should be one of \code{functions}, \code{data}, \code{transformed data}, \code{parameters}, \code{transformed parameters}, \code{model}, or \code{generated quantities}. The named \code{section} slot should be one of \code{declare} or \code{body}, specifying whether the new code appears in the declared variables header or rather the code body of a block. The named \code{pos} slot should be one of \code{top}, \code{bottom}, or \code{pattern}. The position \code{pattern} uses the additional named slot \code{pattern} to search-and-replace, replacing the text in \code{pattern} with the text in the first argument. See the example at the end of this help page.
  
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

