rethinking
==========

This R package accompanies a course and book on Bayesian data analysis. It contains tools for conducting both MAP estimation and Hamiltonian Monte Carlo (through RStan). These tools force the user to specify the model as a list of explicit distributional assumptions.

For example, a simple linear regression could be specified with this list of formulas:

```
f <- alist(
    y ~ dnorm( mu , sigma ),
    mu ~ dnorm( 0 , 10 ),
    sigma ~ dcauchy( 0 , 1 )
)
```

The first formula in the list is the likelihood; the second is the prior for ``mu``; the third is the prior for ``sigma`` (implicitly a half-Cauchy, due to positive constraint on ``sigma``).

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

The same formula list can be compiled into a Stan (mc-stan.org) model:

```
fit.stan <- map2stan( 
    f , 
    data=list(y=c(-1,1)) , 
    start=list(mu=0,sigma=1)
)
```

The chain runs automatically, provided ``rstan`` is installed. The Stan code can be accessed by using ``stancode(fit.stan)``:

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
And ``map2stan`` supports decomposition of covariance matrices into vectors of standard deviations and a correlation matrix, such that priors can be specified independently for each:
```
f4 <- alist(
    y ~ dnorm( mu , sigma ),
    mu ~ a + aj + (b + bj)*x,
    c(aj,bj)[group] <- dmvnorm2( 0 , sigma_group , Rho_group ),
    a ~ dnorm( 0 , 10 ),
    b ~ dnorm( 0 , 1 ),
    sigma ~ dcauchy( 0 , 1 ),
    sigma_group ~ dcauchy( 0 , 1 ),
    Rho_group ~ dlkjcorr(2)
)
```
