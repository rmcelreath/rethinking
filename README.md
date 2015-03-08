rethinking
==========

This R package accompanies a course and book on Bayesian data analysis. It contains tools for conducting both MAP estimation and Hamiltonian Monte Carlo (through RStan). These tools force the user to specify the model as a list of explicit distributional assumptions.

For example, a simple Gaussian model could be specified with this list of formulas:

```
f <- alist(
    y ~ dnorm( mu , sigma ),
    mu ~ dnorm( 0 , 10 ),
    sigma ~ dcauchy( 0 , 1 )
)
```

The first formula in the list is the likelihood; the second is the prior for ``mu``; the third is the prior for ``sigma`` (implicitly a half-Cauchy, due to positive constraint on ``sigma``).

## MAP estimation

Then to use maximum a posteriori (MAP) fitting:
```
library(rethinking)
fit <- map( 
    f , 
    data=list(y=c(-1,1)) , 
    start=list(mu=0,sigma=1)
)
```
The object ``fit`` holds the result.

## Hamiltonian Monte Carlo estimation

The same formula list can be compiled into a Stan (mc-stan.org) model:

```
fit.stan <- map2stan( 
    f , 
    data=list(y=c(-1,1)) , 
    start=list(mu=0,sigma=1)
)
```
The ``start`` list is optional, provided a prior is defined for every parameter. In that case, ``map2stan`` will automatically sample from each prior to get starting values for the chains. The chain runs automatically, provided ``rstan`` is installed. The Stan code can be accessed by using ``stancode(fit.stan)``:
```
data{
    int<lower=1> N;
    real y[N];
}
parameters{
    real mu;
    real<lower=0> sigma;
}
model{
    mu ~ normal( 0 , 10 );
    sigma ~ cauchy( 0 , 1 );
    y ~ normal( mu , sigma );
}
generated quantities{
    real dev;
    dev <- 0;
    dev <- dev + (-2)*normal_log( y , mu , sigma );
}
```

## Posterior prediction

Both ``map`` and ``map2stan`` model fits can be post-processed to produce posterior distributions of any linear models and posterior predictive distributions.

``link`` is used to compute values of any linear models over samples from the posterior distribution. 

``sim`` is used to simulate posterior predictive distributions, simulating outcomes over samples from the posterior distribution of parameters. See ``?link`` and ``?sim`` for details.

``postcheck`` automatically computes posterior predictive (retrodictive?) checks for each case used to fit a model. 


## Multilevel model formulas

While ``map`` is limited to fixed effects models for the most part, ``map2stan`` can specify multilevel models, even quite complex ones. For example, a simple varying intercepts model looks like:
```
f2 <- alist(
    y ~ dnorm( mu , sigma ),
    mu <- a + aj,
    aj[group] ~ dnorm( 0 , sigma_group ),
    a ~ dnorm( 0 , 10 ),
    sigma ~ dcauchy( 0 , 1 ),
    sigma_group ~ dcauchy( 0 , 1 )
)
```
And with varying slopes as well:
```
f3 <- alist(
    y ~ dnorm( mu , sigma ),
    mu <- a + aj + (b + bj)*x,
    c(aj,bj)[group] ~ dmvnorm( 0 , Sigma_group ),
    a ~ dnorm( 0 , 10 ),
    b ~ dnorm( 0 , 1 ),
    sigma ~ dcauchy( 0 , 1 ),
    Sigma_group ~ inv_wishart( 3 , diag(2) )
)
```

## Nice covariance priors

And ``map2stan`` supports decomposition of covariance matrices into vectors of standard deviations and a correlation matrix, such that priors can be specified independently for each:
```
f4 <- alist(
    y ~ dnorm( mu , sigma ),
    mu <- a + aj + (b + bj)*x,
    c(aj,bj)[group] ~ dmvnorm2( 0 , sigma_group , Rho_group ),
    a ~ dnorm( 0 , 10 ),
    b ~ dnorm( 0 , 1 ),
    sigma ~ dcauchy( 0 , 1 ),
    sigma_group ~ dcauchy( 0 , 1 ),
    Rho_group ~ dlkjcorr(2)
)
```

## Semi-automated Bayesian imputation

It is possible to code simple Bayesian imputations this way. For example, let's simulate a simple regression with missing predictor values:
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
What ``map2stan`` does is notice the missing values, see the distribution assigned to the variable with the missing values, build the Stan code that uses a mix of observed and estimated ``x`` values in the regression. See the ``stancode(m)`` for details of the implementation.

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
d$island <- 1:10
mGP <- map2stan(
    alist(
        total_tools ~ dpois( mu ),
        log(mu) <- a + aj[island],
        a ~ dnorm(0,10),
        aj[island] ~ GPL2( Dmat , etasq , rhosq , 0.01 ),
        etasq ~ dcauchy(0,1),
        rhosq ~ dcauchy(0,1)
    ),
    data=list(
        total_tools=d$total_tools,
        island=d$island,
        Dmat=islandsDistMatrix),
    constraints=list(
        etasq="lower=0",
        rhosq="lower=0"
    ),
    warmup=1000 , iter=5000 , chains=4 )
```
Note the use of the ``constraints`` list to pass custom parameter constraints to Stan. This example is explored in more detail in the (in prep) book. 

## Information criteria

Both ``map`` and ``map2stan`` provide DIC and WAIC. Well, in most cases they do. In truth, both tools are flexible enough that you can specify models for which neither DIC nor WAIC can be correctly calculated. But for ordinary GLMs and GLMMs, it works. See the R help ``?WAIC``. A convenience function ``compare`` summarizes information criteria comparisons, including standard errors for WAIC. 

``ensemble`` computes ``link`` and ``sim`` output for an ensemble of models, each weighted by its Akaike weight, as computed from WAIC.

