rethinking
==========

This R package accompanies a course and book on Bayesian data analysis: McElreath 2020. Statistical Rethinking, 2nd edition, CRC Press. If you are using it with the first edition of the book, please see the notes at the bottom of this file.

It contains tools for conducting both quick quadratic approximation of the posterior distribution as well as Hamiltonian Monte Carlo (through RStan or cmdstanr - mc-stan.org). Many packages do this. The signature difference of this package is that it forces the user to specify the model as a list of explicit distributional assumptions. This is more tedious than typical formula-based tools, but it is also much more flexible and powerful and---most important---useful for teaching and learning. When students have to write out every detail of the model, they actually learn the model.

For example, a simple Gaussian model could be specified with this list of formulas:

```
f <- alist(
    y ~ dnorm( mu , sigma ),
    mu ~ dnorm( 0 , 10 ),
    sigma ~ dexp( 1 )
)
```

The first formula in the list is the probability of the outcome (likelihood); the second is the prior for ``mu``; the third is the prior for ``sigma``.

# Installation

There are three steps. (1) Install ``rstan``, (2) install ``cmdstanr``, (3) install ``rethinking``. Details follow.

First, install the C++ toolchain and install the ``rstan`` package. Go to ``https://mc-stan.org/users/interfaces/rstan.html`` and follow the instructions for your platform. The biggest challenge is getting a C++ compiler configured to work with your installation of R. The instructions are quite thorough. Obey them, and you'll succeed.

Second, install the ``cmdstanr`` package. Visit ``https://mc-stan.org/cmdstanr/``. The first time you install cmdstanr, you will also need compile the libraries with ``cmdstanr::install_cmdstan()``. All this of this bother is worth it. You just have to do it once.

Third, once rstan and cmdstanr are installed (almost there), then you can install ``rethinking`` from within R using:
```
install.packages(c("coda","mvtnorm","devtools","loo","dagitty","shape"))
devtools::install_github("rmcelreath/rethinking")
```
If there are any problems, they likely arise when trying to install ``rstan``, so the ``rethinking`` package has little to do with it. See the manual linked above for some hints about getting ``rstan`` installed. But always consult the RStan section of the website at ``mc-stan.org`` for the latest information on RStan.

Note that the ``rethinking`` package is not on CRAN, just on github. The ``rethinking`` package is never going to be on CRAN.

# rethinking slim - no MCMC

If you just want to work through the first half of the course, without bothering with MCMC and Stan installs, you can install the 'slim' version of the rethinking package. Do this:
```
install.packages(c("coda","mvtnorm","devtools","loo","dagitty"))
devtools::install_github("rmcelreath/rethinking@slim")
```
The ``quap`` function and related helper functions should still work, and you'll be able to work through Chapter 8 before you need to install the full version with Stan.

# Quadratic Approximation with `quap`

Almost any ordinary generalized linear model can be specified with `quap`. To use quadratic approximation:
```
library(rethinking)

f <- alist(
    y ~ dnorm( mu , sigma ),
    mu ~ dnorm( 0 , 10 ),
    sigma ~ dexp( 1 )
)

fit <- quap( 
    f , 
    data=list(y=c(-1,1)) , 
    start=list(mu=0,sigma=1)
)
```
The object ``fit`` holds the result. For a summary of marginal posterior distributions, use ``summary(fit)`` or ``precis(fit)``:
```
      mean   sd  5.5% 94.5%
mu    0.00 0.59 -0.95  0.95
sigma 0.84 0.33  0.31  1.36
```
It also supports vectorized parameters, which is convenient for categories. See examples `?quap`. 

In the first edition of the textbook, this function was called `map`. It can still be used with that alias. It was renamed, because the name `map` was misleading. This function produces quadratic approximations of the posterior distribution, not just maximum a posteriori (MAP) estimates.

# Hamiltonian Monte Carlo with `ulam` (and `map2stan`)

The same formula list can be compiled into a Stan (mc-stan.org) model using one of two tools: `ulam` or `map2stan`. For simple models, they are identical. `ulam` is the newer tool that allows for much more flexibility, including explicit variable types and custom distributions. `map2stan` is the original tool from the first edition of the package and textbook. Going forward, new features will be added to `ulam`.

`ulam` is named after Stanisław Ulam, who was one of the parents of the Monte Carlo method and is the namesake of the Stan project as well. It is pronounced something like [OO-lahm], not like [YOU-lamm].

Both tools take the same kind of input as `quap`:
```
fit_stan <- ulam( f , data=list(y=c(-1,1)) )
```
The chain runs automatically, provided ``rstan`` is installed. Chain diagnostics are displayed in the `precis(fit_stan)` output:
```
      mean   sd  5.5% 94.5% n_eff Rhat
sigma 1.45 0.72  0.67  2.84   145    1
mu    0.12 1.04 -1.46  1.59   163    1
```
For `ulam` models, `plot` displays the same information as `precis` and `traceplot` displays the chains.

`extract.samples` returns samples in a list. `extract.prior` samples from the prior and returns the samples in a list as well.

The `stanfit` object itself is in the `@stanfit` slot. Anything you'd do with a Stan model can be done with that slot directly.

The Stan code can be accessed by using ``stancode(fit_stan)``:
```
data{
    real y[2];
}
parameters{
    real<lower=0> sigma;
    real mu;
}
model{
    sigma ~ exponential( 1 );
    mu ~ normal( 0 , 10 );
    y ~ normal( mu , sigma );
}
```

Note that `ulam` doesn't care about R distribution names. You can instead use Stan-style names:
```
fit_stan <- ulam(
    alist(
        y ~ normal( mu , sigma ),
        mu ~ normal( 0 , 10 ),
        sigma ~ exponential( 1 )
    ), data=list(y=c(-1,1)) )
```

## Posterior prediction

All ``quap``, `ulam`, and ``map2stan`` objects can be post-processed to produce posterior predictive distributions.

``link`` is used to compute values of any linear models over samples from the posterior distribution. 

``sim`` is used to simulate posterior predictive distributions, simulating outcomes over samples from the posterior distribution of parameters. `sim` can also be used to simulate prior predictives.

See ``?link`` and ``?sim`` for details.

``postcheck`` automatically computes posterior predictive (retrodictive?) checks. It merely uses `link` and `sim`.


## Multilevel model formulas

While ``quap`` is limited to fixed effects models for the most part, `ulam` can specify multilevel models, even quite complex ones. For example, a simple varying intercepts model looks like:
```
# prep data
data( UCBadmit )
UCBadmit$male <- as.integer(UCBadmit$applicant.gender=="male")
UCBadmit$dept <- rep( 1:6 , each=2 )
UCBadmit$applicant.gender <- NULL

# varying intercepts model
m_glmm1 <- ulam(
    alist(
        admit ~ binomial(applications,p),
        logit(p) <- a[dept] + b*male,
        a[dept] ~ normal( abar , sigma ),
        abar ~ normal( 0 , 4 ),
        sigma ~ half_normal(0,1),
        b ~ normal(0,1)
    ), data=UCBadmit )
```
The analogous varying slopes model is:
```
m_glmm2 <- ulam(
    alist(
        admit ~ binomial(applications,p),
        logit(p) <- a[dept] + b[dept]*male,
        c( a , b )[dept] ~ multi_normal( c(abar,bbar) , Rho , sigma ),
        abar ~ normal( 0 , 4 ),
        bbar ~ normal(0,1),
        sigma ~ half_normal(0,1),
        Rho ~ lkjcorr(2)
    ),
    data=UCBadmit )
```
Another way to express the varying slopes model is with a vector of varying effects. This is made possible by using an explicit vector declaration inside the formula:
```
m_glmm3 <- ulam(
    alist(
        admit ~ binomial(applications,p),
        logit(p) <- v[dept,1] + v[dept,2]*male,
        vector[2]:v[dept] ~ multi_normal( c(abar,bbar) , Rho , sigma ),
        abar ~ normal( 0 , 4 ),
        bbar ~ normal(0,1),
        sigma ~ half_normal(0,1),
        Rho ~ lkjcorr(2)
    ),
    data=UCBadmit )
```
That `vector[2]:v[dept]` means "declare a vector of length two for each unique dept". To access the elements of these vectors, the linear model uses multiple indexes inside the brackets: `[dept,1]`. 

This strategy can be taken one step further and the means can be declared as a vector as well:
```
m_glmm4 <- ulam(
    alist(
        admit ~ binomial(applications,p),
        logit(p) <- v[dept,1] + v[dept,2]*male,
        vector[2]:v[dept] ~ multi_normal( v_mu , Rho , sigma ),
        vector[2]:v_mu ~ normal(0,1),
        sigma[1] ~ half_normal(0,1),
        sigma[2] ~ half_normal(0,2),
        Rho ~ lkjcorr(2)
    ),
    data=UCBadmit )
```
And a completely non-centered parameterization can be coded directly as well:
```
m_glmm5 <- ulam(
    alist(
        admit ~ binomial(applications,p),
        logit(p) <- v_mu[1] + v[dept,1] + (v_mu[2] + v[dept,2])*male,
        matrix[dept,2]: v <- t(diag_pre_multiply( sigma , L_Rho ) * z),
        matrix[2,dept]: z ~ normal( 0 , 1 ),
        vector[2]: v_mu[[1]] ~ normal(0,4),
        vector[2]: v_mu[[2]] ~ normal(0,1),
        vector[2]: sigma ~ half_normal(0,1),
        cholesky_factor_corr[2]: L_Rho ~ lkj_corr_cholesky( 2 )
    ),
    data=UCBadmit )
```
In the above, the varying effects matrix `v` is constructed from a matrix of z-scores `z` and a covariance structure contained in `sigma` and a Cholesky factor `L_Rho`. Note the double-bracket notation `v_mu[[1]]` allowing distinct priors for each index of a vector.

## log-likelihood calculations for WAIC and LOOCV

`ulam` can optionally return pointwise log-likelihood values. These are needed for computing WAIC and PSIS-LOO. The `log_lik` argument toggles this on:
```
m_glmm1 <- ulam(
    alist(
        admit ~ binomial(applications,p),
        logit(p) <- a[dept] + b*male,
        a[dept] ~ normal( abar , sigma ),
        abar ~ normal( 0 , 4 ),
        sigma ~ half_normal(0,1),
        b ~ normal(0,1)
    ), data=UCBadmit , log_lik=TRUE )
WAIC(m_glmm1)
```
The additional code has been added to the generated quantities block of the Stan model (see this with `stancode(m_glmm1)`):
```
generated quantities{
    vector[12] log_lik;
    vector[12] p;
    for ( i in 1:12 ) {
        p[i] = a[dept[i]] + b * male[i];
        p[i] = inv_logit(p[i]);
    }
    for ( i in 1:12 ) log_lik[i] = binomial_lpmf( admit[i] | applications[i] , p[i] );
}
```

## Conditional statements, custom distributions, and mixture models

`ulam` also supports if-then statements and custom distribution assignments. These are useful for coding mixture models, such as zero-inflated Poisson and discrete missing value models.

Here's an example zero-inflated Poisson model.
```
# zero-inflated poisson
# gen data first - example from text
prob_drink <- 0.2 # 20% of days
rate_work <- 1    # average 1 manuscript per day
N <- 365
drink <- rbinom( N , 1 , prob_drink )
y <- as.integer( (1-drink)*rpois( N , rate_work ) )
x <- rnorm( N ) # dummy covariate

# now ulam code
m_zip <- ulam(
    alist(
        y|y==0 ~ custom( log_mix( p , 0 , poisson_lpmf(0|lambda) ) ),
        y|y>0 ~ custom( log1m(p) + poisson_lpmf(y|lambda) ),
        logit(p) <- ap,
        log(lambda) <- al + bl*x,
        ap ~ dnorm(0,1),
        al ~ dnorm(0,10),
        bl ~ normal(0,1)
    ) ,
    data=list(y=y,x=x) )
```
The Stan code corresponding to the first two lines in the formula above is:
```
for ( i in 1:365 ) 
    if ( y[i] > 0 ) target += log1m(p) + poisson_lpmf(y[i] | lambda[i]);
for ( i in 1:365 ) 
    if ( y[i] == 0 ) target += log_mix(p, 0, poisson_lpmf(0 | lambda[i]));
```
What `custom` does is define custom `target` updates. And the `|` operator makes the line conditional. Note that `log1m`, `log_mix`, and `poisson_lpmf` are Stan functions.

The same `custom` distribution approach allows for marginalization over discrete missing values. Let's introduce some missing values in the `UCBadmit` data from earlier.
```
UCBadmit$male2 <- UCBadmit$male
UCBadmit$male2[1:2] <- (-1) # missingness code
UCBadmit$male2 <- as.integer(UCBadmit$male2)
```
Now the model needs to detect when `male2` is missing (-1) and then compute a mixture over the unknown state.
```
m_mix <- ulam(
    alist(
        admit|male2==-1 ~ custom( log_mix( 
            phi_male , 
            binomial_lpmf(admit|applications,p_m1) , 
            binomial_lpmf(admit|applications,p_m0) ) ),
        admit|male2>-1 ~ binomial( applications , p ),
        logit(p) <- a[dept] + b*male2,
        logit(p_m1) <- a[dept] + b*1,
        logit(p_m0) <- a[dept] + b*0,
        male2|male2>-1 ~ bernoulli( phi_male ),
        phi_male ~ beta(2,2),
        a[dept] ~ normal(0,4),
        b ~ normal(0,1)
    ),
    data=UCBadmit )
```
Note the addition of `phi_male` to average over the unknown state. 

## Continuous missing data imputation

In principle, imputation of missing real-valued data is easy: Just replace each missing value with a parameter. In practice, this involves a bunch of annoying bookkeeping. `ulam` has a macro named `merge_missing` to simplify this.
```
UCBadmit$x <- rnorm(12)
UCBadmit$x[1:2] <- NA
m_miss <- ulam(
    alist(
        admit ~ binomial(applications,p),
        logit(p) <- a + b*male + bx*x_merge,
        x_merge ~ normal( 0 , 1 ),
        x_merge <- merge_missing( x , x_impute ),
        a ~ normal(0,4),
        b ~ normal(0,1),
        bx ~ normal(0,1)
    ),
    data=UCBadmit )
```
What `merge_missing` does is find the `NA` values in `x` (whichever symbol is the first argument), build a vector of parameters called `x_impute` (whatever you name the second argument) of the right length, and piece together a vector `x_merge` that contains both, in the right places. You can then assign a prior to this vector and use it in linear models as usual.

The merging is done as the Stan model runs, using a custom function block. See the Stan code `stancode(m_miss)` for all the lovely details.

`merge missing` is an example of a macro, which is a way for `ulam` to use function names to trigger special compilation. In this case, `merge_missing` both inserts a function in the Stan model and builds the necessary index to locate the missing values during run time. Macros will get full documentation later, once the system is finalized.

## Gaussian processes

A simple Gaussian process, like the Oceanic islands example in Chapter 13 of the book, is done as:
```
data(Kline2)
d <- Kline2
data(islandsDistMatrix)
d$society <- 1:10
dat <- list(
    y=d$total_tools,
    society=d$society,
    log_pop = log(d$population),
    Dmat=islandsDistMatrix
)

m_GP1 <- ulam(
    alist(
        y ~ poisson( mu ),
        log(mu) <- a + aj[society] + b*log_pop,
        a ~ normal(0,10),
        b ~ normal(0,1),
        vector[10]: aj ~ multi_normal( 0 , SIGMA ),
        matrix[10,10]: SIGMA <- cov_GPL2( Dmat , etasq , rhosq , 0.01 ),
        etasq ~ exponential(1),
        rhosq ~ exponential(1)
    ),
    data=dat )
```
This is just an ordinary varying intercepts model, but all 10 intercepts are drawn from a single Gaussian distribution. The covariance matrix `SIGMA` is defined in the usual L2-norm. Again, `cov_GPL2` is a macro that inserts a function in the Stan code to compute the covariance matrix as the model runs.

Fancier Gaussian processes require a different parameterization. And these can be built as well. Here's an example using 151 primate species and a phylogenetic distance matrix. First, prepare the data:
```
data(Primates301)
data(Primates301_distance_matrix)
d <- Primates301
d$name <- as.character(d$name)
dstan <- d[ complete.cases( d$social_learning, d$research_effort , d$body , d$brain ) , ]
# prune distance matrix to spp in dstan
spp_obs <- dstan$name
y <- Primates301_distance_matrix
y2 <- y[ spp_obs , spp_obs ]
# scale distances
y3 <- y2/max(y2)
```
Now the model, which is a non-centered L2-norm Gaussian process:
```
m_GP2 <- ulam(
    alist(
        social_learning ~ poisson( lambda ),
        log(lambda) <- a + g[spp_id] + b_ef*log_research_effort + b_body*log_body + b_eq*log_brain,
        a ~ normal(0,1),
        vector[N_spp]: g <<- L_SIGMA * eta,
        vector[N_spp]: eta ~ normal( 0 , 1 ),
        matrix[N_spp,N_spp]: L_SIGMA <<- cholesky_decompose( SIGMA ),
        matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat , etasq , rhosq , 0.01 ),
        b_body ~ normal(0,1),
        b_eq ~ normal(0,1),
        b_ef ~ normal(1,1),
        etasq ~ exponential(1),
        rhosq ~ exponential(1)
    ),
    data=list(
        N_spp = nrow(dstan),
        social_learning = dstan$social_learning,
        spp_id = 1:nrow(dstan),
        log_research_effort = log(dstan$research_effort),
        log_body = log(dstan$body),
        log_brain = log(dstan$brain),
        Dmat = y3
    ) , 
    control=list(max_treedepth=15,adapt_delta=0.95) ,
    sample=FALSE )
```
This model does not sample quickly, so I've set `sample=FALSE`. You can still inspect the Stan code with `stancode(m_GP2)`. 

Note that the covariance `SIGMA` is built the same way as before, but then we immediately decompose it to a Cholesky factor and build the varying intercepts `g` by matrix multiplication. The `<<-` operator tells `ulam` not to loop, but to do a direct assignment. So `g <<- L_SIGMA * eta` does the right linear algebra.

## Within-chain multithreading

Using ``cmdstanr`` instead of ``rstan`` is currently the only way to use within-chain multithreading with ``rethinking``. It also tends to compile models faster and is more intelligent about when models need to be re-compiled, so using ``cmdstanr`` is recommended, even if you don't want multithreading.

If you want ``ulam`` to access Stan using the ``cmdstanr`` package, then you may install that as well with
```
devtools::install_github("stan-dev/cmdstanr")
```
If you haven't installed cmdstan previously, you will also need to do that with ``install_cmdstan()``. 

Then you need to add ``cmdstan=TRUE`` to the ``ulam`` code. The ``threads`` argument controls the number of threads per chain. Example:
```
N <- 1e4
x <- rnorm(N)
m <- 1 + rpois(N,2)
y <- rbinom( N , size=m , prob=inv_logit(-3+x) )
dat <- list( y=y , x=x , m=m )
# two threads
m1 <- ulam(
    alist(
        y ~ binomial_logit( m , logit_p ),
        logit_p <- a + b*x,
        a ~ normal(0,1.5),
        b ~ normal(0,0.5)
    ) , data=dat , 
    cmdstan=TRUE , threads=2 , refresh=1000 )
```
There are models that cannot be automaticaly multithreaded this way, because of the complexity of the code. In those cases, you can write the code directly in Stan. See [this guide](https://mc-stan.org/users/documentation/case-studies/reduce_sum_tutorial.html). Writing multithreaded models direct in Stan can also be more efficient, since you can make detailed choices about which variables to pass and which pieces of the model to multithread.

## Work in progress

`ulam` is still in development, but mostly feature complete. It will remain primarily a teaching tool, exposing the statistical details of the model while hiding some of the programming details necessary in Stan.


# `map2stan` syntax and features

The older `map2stan` function makes stronger assumtions about the formulas it will see. This allows is to provide some additional automation and it has some special syntax as a result. `ulam` in contrast supports such features through its macros library.

## Non-centered parameterization
Here is a non-centered parameterization that moves the scale parameters in the varying effects prior to the linear model, which is often more efficient for sampling:
```
f4u <- alist(
    y ~ dnorm( mu , sigma ),
    mu <- a + zaj[group]*sigma_group[1] + 
         (b + zbj[group]*sigma_group[2])*x,
    c(zaj,zbj)[group] ~ dmvnorm( 0 , Rho_group ),
    a ~ dnorm( 0 , 10 ),
    b ~ dnorm( 0 , 1 ),
    sigma ~ dcauchy( 0 , 1 ),
    sigma_group ~ dcauchy( 0 , 1 ),
    Rho_group ~ dlkjcorr(2)
)
```
Chapter 13 of the book provides a lot more detail on this issue. 

We can take this strategy one step further and remove the correlation matrix, ``Rho_group``, from the prior as well.  ``map2stan`` facilitates this form via the ``dmvnormNC`` density, which uses an internal Cholesky decomposition of the correlation matrix to build the varying effects. Here is the previous varying slopes model, now with the non-centered notation:
```
f4nc <- alist(
    y ~ dnorm( mu , sigma ),
    mu <- a + aj[group] + (b + bj[group])*x,
    c(aj,bj)[group] ~ dmvnormNC( sigma_group , Rho_group ),
    a ~ dnorm( 0 , 10 ),
    b ~ dnorm( 0 , 1 ),
    sigma ~ dcauchy( 0 , 1 ),
    sigma_group ~ dcauchy( 0 , 1 ),
    Rho_group ~ dlkjcorr(2)
)
```
Internally, a Cholesky factor ``L_Rho_group`` is used to perform sampling. It will appear in the returned samples, in addition to ``Rho_group``, which is constructed from it.

## Semi-automated Bayesian imputation

It is possible to code simple Bayesian imputations. For example, let's simulate a simple regression with missing predictor values:
```
N <- 100
N_miss <- 10
x <- rnorm( N )
y <- rnorm( N , 2*x , 1 )
x[ sample(1:N,size=N_miss) ] <- NA
```
That removes 10 ``x`` values. Then the ``map2stan`` formula list just defines a distribution for ``x``:
```
f5 <- alist(
    y ~ dnorm( mu , sigma ),
    mu <- a + b*x,
    x ~ dnorm( mu_x, sigma_x ),
    a ~ dnorm( 0 , 100 ),
    b ~ dnorm( 0  , 10 ),
    mu_x ~ dnorm( 0 , 100 ),
    sigma_x ~ dcauchy(0,2),
    sigma ~ dcauchy(0,2)
)
m5 <- map2stan( f5 , data=list(y=y,x=x) )
```
What ``map2stan`` does is notice the missing values, see the distribution assigned to the variable with the missing values, build the Stan code that uses a mix of observed and estimated ``x`` values in the regression. See the ``stancode(m5)`` for details of the implementation.

## Semi-automated marginalization for binary discrete missing values

Binary (0/1) variables with missing values present a special obstacle, because Stan cannot sample discrete parameters. So instead of imputing binary missing values, ``map2stan`` can average (marginalize) over them. As in the above case, when ``map2stan`` detects missing values in a predictor variable, it will try to find a distribution for the variable containing them. If this variable is binary (0/1), then it will construct a mixture model in which each term is the log-likelihood conditional on the variables taking a particular combination of 0/1 values.

Following the example in the previous section, we can simulate missingness in a binary predictor:
```
N <- 100
N_miss <- 10
x <- rbinom( N , size=1 , prob=0.5 )
y <- rnorm( N , 2*x , 1 )
x[ sample(1:N,size=N_miss) ] <- NA
```
The model definition is analogous to the previous, but also requires some care in specifying constraints for the hyperparameters that define the distribution for ``x``:
```
f6 <- alist(
    y ~ dnorm( mu , sigma ),
    mu <- a + b*x,
    x ~ bernoulli( phi ),
    a ~ dnorm( 0 , 100 ),
    b ~ dnorm( 0  , 10 ),
    phi ~ beta( 1 , 1 ),
    sigma ~ dcauchy(0,2)
)
m6 <- map2stan( f6 , data=list(y=y,x=x) , constraints=list(phi="lower=0,upper=1") )
```
The algorithm works, in theory, for any number of binary predictors with missing values. For example, with two predictors, each with missingness:
```
N <- 100
N_miss <- 10
x1 <- rbinom( N , size=1 , prob=0.5 )
x2 <- rbinom( N , size=1 , prob=0.1 )
y <- rnorm( N , 2*x1 - x2  , 1 )
x1[ sample(1:N,size=N_miss) ] <- NA
x2[ sample(1:N,size=N_miss) ] <- NA
f7 <- alist(
    y ~ dnorm( mu , sigma ),
    mu <- a + b1*x1 + b2*x2,
    x1 ~ bernoulli( phi1 ),
    x2 ~ bernoulli( phi2 ),
    a ~ dnorm( 0 , 100 ),
    c(b1,b2) ~ dnorm( 0  , 10 ),
    phi1 ~ beta( 1 , 1 ),
    phi2 ~ beta( 1 , 1 ),
    sigma ~ dcauchy(0,2)
)
m7 <- map2stan( f7 , data=list(y=y,x1=x1,x2=x2) , 
      constraints=list(phi1="lower=0,upper=1",phi2="lower=0,upper=1") )
```
While the unobserved values for the binary predictors are usually not of interest, they can be computed from the posterior distribution. Adding the argument ``do_discrete_imputation=TRUE`` instructs ``map2stan`` to perform these calculations automatically. Example:
```
m6 <- map2stan( f6 , data=list(y=y,x=x) , constraints=list(phi="lower=0,upper=1") ,
      do_discrete_imputation=TRUE )
precis( m6 , depth=2 )
```
The output contains samples for each case with imputed probilities that ``x`` takes the value 1.

The algorithm works by constructing a list of mixture terms that are needed to to compute the probability of each observed ``y`` value. In the simplest case, with only one predictor with missing values, the implied mixture likelihood contains two terms:
```
Pr(y[i]) = Pr(x[i]=1)Pr(y[i]|x[i]=1) + Pr(x[i]=0)Pr(y[i]|x[i]=0)
```
In the parameters of our example model ``m6`` above, this is:
```
Pr(y[i]) = phi*N(y[i]|a+b,sigma) + (1-phi)*N(y[i]|a,sigma)
```
It is now a simple matter to loop over cases ``i`` and compute the above for each. Similarly the posterior probability of that ``x[i]==1`` is given as:
```
Pr(x[i]==1|y[i]) = phi*N(y[i]|a+b,sigma) / Pr(y[i])
```
When only one predictor has missingness, then this is simple. What about when there are two or more? In that case, all the possible combinations of missingness have to be accounted for. For example, suppose there are two predictors, ``x1`` and ``x2``, both with missingness on case ``i``. Now the implied mixture likelihood is:
```
Pr(y[i]) = Pr(x1=1)Pr(x2=1)*Pr(y[i]|x1=1,x2=1) + Pr(x1=1)Pr(x2=0)Pr(y[i]|x1=1,x2=0) + Pr(x1=0)Pr(x2=1)Pr(y[i]|x1=0,x2=1) + Pr(x1=0)Pr(x2=0)Pr(y[i]|x1=0,x2=0)
```
There are four combinations of unobserved values, and so four terms in the mixture likelihood. When ``x2`` is instead observed, we can substitute the observed value into the above, and then the mixture simplifies readily to our previous two-term likelihood:
```
Pr(y[i]|x2[i]==1) = Pr(x1=1)Pr(x2=1)Pr(y[i]|x1=1,x2=1) + Pr(x1=1)Pr(x2=0)Pr(y[i]|x1=1,x2=1) + Pr(x1=0)Pr(x2=1)Pr(y[i]|x1=0,x2=1) + Pr(x1=0)Pr(x2=0)Pr(y[i]|x1=0,x2=1)
                  = [Pr(x1=1)Pr(x2=1)+Pr(x1=1)Pr(x2=0)]Pr(y[i]|x1=1,x2=1) 
                    + [Pr(x1=0)Pr(x2=1)+Pr(x1=0)Pr(x2=0)]Pr(y[i]|x1=0,x2=1)
                  = Pr(x1=1)Pr(y[i]|x1=1,x2=1) + Pr(x1=0)Pr(y[i]|x1=0,x2=1)
```
This implies that if we loop over cases ``i`` and insert any observed values into the general mixture likelihood, we can compute the relevant mixture for the specific combination of missingness on each case ``i``. That is what ``map2stan`` does. The general mixture terms can be generated algorithmically. The code below generates a matrix of terms for ``n`` binary variables with missingness.
```
ncombinations <- 2^n
d <- matrix(NA,nrow=ncombinations,ncol=n)
for ( col_var in 1:n ) 
    d[,col_var] <- rep( 0:1 , each=2^(col_var-1) , length.out=ncombinations )
```
Rows of ``d`` contain terms, columns contain variables, and the values in each column are the corresponding values of each variable. The algorithm builds a linear model for each row in this matrix, composes the mixture likelihood as the sum of these rows, and performs proper substitutions of observed values. All calculations are done on the log scale, for precision.

## Gaussian process

A basic Gaussian process can be specified with the ``GPL2`` distribution label. This implies a multivariate Gaussian with a covariance matrix defined by the ordinary L2 norm distance function:
```
k(i,j) = eta^2 * exp( -rho^2 * D(i,j)^2 ) + ifelse(i==j,sigma^2,0)
```
where ``D`` is a matrix of pairwise distances. To use this convention in, for example, a spatial autocorrelation model:
```
library(rethinking)
data(Kline2)
d <- Kline2
data(islandsDistMatrix)
d$society <- 1:10
mGP <- map2stan(
    alist(
        total_tools ~ dpois( mu ),
        log(mu) <- a + aj[society],
        a ~ dnorm(0,10),
        aj[society] ~ GPL2( Dmat , etasq , rhosq , 0.01 ),
        etasq ~ dcauchy(0,1),
        rhosq ~ dcauchy(0,1)
    ),
    data=list(
        total_tools=d$total_tools,
        society=d$society,
        Dmat=islandsDistMatrix),
    constraints=list(
        etasq="lower=0",
        rhosq="lower=0"
    ),
    warmup=1000 , iter=5000 , chains=4 )
```
Note the use of the ``constraints`` list to pass custom parameter constraints to Stan. This example is explored in more detail in the book. 

## Information criteria

Both ``map`` and ``map2stan`` provide DIC and WAIC. Well, in most cases they do. In truth, both tools are flexible enough that you can specify models for which neither DIC nor WAIC can be correctly calculated. But for ordinary GLMs and GLMMs, it works. See the R help ``?WAIC``. A convenience function ``compare`` summarizes information criteria comparisons, including standard errors for WAIC. 

`ulam` supports WAIC calculation with the optional `log_lik=TRUE` argument, which returns the kind of log-likelihood vector needed by the `loo` package.

``ensemble`` computes ``link`` and ``sim`` output for an ensemble of models, each weighted by its Akaike weight, as computed from WAIC.

# Code issues with 1st edition of Statistical Rethinking

A small change to ``link`` has broken two examples in the first edition of the book, in Chapter 7.

## R code 7.10
> mu.Africa.mean <- apply( mu.Africa , 2 , mean )
Error in apply(mu.Africa, 2, mean) : dim(X) must have a positive length

This occurs because link() now returns all linear models. So mu.Africa is a list containing mu and gamma. To fix, use:

> mu.Africa.mean <- apply( mu.Africa$mu , 2 , mean )

Use a similar fix in the other apply() calls in the same section.

## R code 7.17
Similar problem as for R code 7.10. Use mu.ruggedlo$mu in place of mu.ruggedlo.

