rethinking
==========

This R package accompanies a course and book on Bayesian data analysis (McElreath 2016. Statistical Rethinking. CRC Press.). It contains tools for conducting both MAP estimation and Hamiltonian Monte Carlo (through RStan - mc-stan.org). These tools force the user to specify the model as a list of explicit distributional assumptions. This is more tedious than typical formula-based tools, but it is also much more flexible and powerful.

For example, a simple Gaussian model could be specified with this list of formulas:

```
f <- alist(
    y ~ dnorm( mu , sigma ),
    mu ~ dnorm( 0 , 10 ),
    sigma ~ dcauchy( 0 , 1 )
)
```

The first formula in the list is the likelihood; the second is the prior for ``mu``; the third is the prior for ``sigma`` (implicitly a half-Cauchy, due to positive constraint on ``sigma``).

## Quick Installation

You can find a manual with expanded installation and usage instructions here: ``http://xcelab.net/rm/software/``

Here's the brief verison. 

You'll need to install ``rstan`` first. Go to ``http://mc-stan.org`` and follow the instructions for your platform. The biggest challenge is getting a C++ compiler configured to work with your installation of R. The instructions at ``https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started`` are quite thorough. Obey them, and you'll likely succeed.

Then you can install ``rethinking`` from within R using:
```
install.packages(c("coda","mvtnorm","devtools","loo"))
library(devtools)
devtools::install_github("rmcelreath/rethinking")
```
If there are any problems, they likely arise when trying to install ``rstan``, so the ``rethinking`` package has little to do with it. See the manual linked above for some hints about getting ``rstan`` installed. But always consult the RStan section of the website at ``mc-stan.org`` for the latest information on RStan.

## MAP estimation

To use maximum a posteriori (MAP) fitting:
```
library(rethinking)

f <- alist(
    y ~ dnorm( mu , sigma ),
    mu ~ dnorm( 0 , 10 ),
    sigma ~ dcauchy( 0 , 1 )
)

fit <- map( 
    f , 
    data=list(y=c(-1,1)) , 
    start=list(mu=0,sigma=1)
)
```
The object ``fit`` holds the result. For a summary of marginal posterior distributions, use ``summary(fit)`` or ``precis(fit)``:
```
      Mean StdDev  2.5% 97.5%
mu    0.00   0.59 -1.16  1.16
sigma 0.84   0.33  0.20  1.48
```

## Hamiltonian Monte Carlo estimation

The same formula list can be compiled into a Stan (mc-stan.org) model:

```
fit.stan <- map2stan( 
    f , 
    data=list(y=c(-1,1)) , 
    start=list(mu=0,sigma=1)
)
```
The ``start`` list is optional, provided a prior is defined for every parameter. In that case, ``map2stan`` will automatically sample from each prior to get starting values for the chains.  The chain runs automatically, provided ``rstan`` is installed. The ``plot`` method will display trace plots for the chains.

The Stan code can be accessed by using ``stancode(fit.stan)``:
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

To run multiple chains in parallel on multiple cores, use the ``cores`` argument:
```
fit.stan <- map2stan( 
    f , 
    data=list(y=c(-1,1)) , 
    start=list(mu=0,sigma=1) ,
    chains=4 , cores=4 , iter=2000 , warmup=1000
)
```
The ``parallel`` package is used here, relying upon ``mclapply`` (Mac, UNIX) or ``parLapply`` (Windows). It is best to run parallel operations in the Terminal/Command Prompt, as GUI interfaces sometimes crash when forking processes.

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
    mu <- a + aj[group],
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
    mu <- a + aj[group] + (b + bj[group])*x,
    c(aj,bj)[group] ~ dmvnorm( 0 , Sigma_group ),
    a ~ dnorm( 0 , 10 ),
    b ~ dnorm( 0 , 1 ),
    sigma ~ dcauchy( 0 , 1 ),
    Sigma_group ~ inv_wishart( 3 , diag(2) )
)
```



## Nice covariance priors

The ``inv_wishart`` prior in the model just above is conventional, but not appealing. Since Stan does not use Gibbs sampling, there is no advantage to the ``inv_wishart`` prior. 
To escape these conventional priors, ``map2stan`` supports decomposition of covariance matrices into vectors of standard deviations and a correlation matrix, such that priors can be specified independently for each:
```
f4 <- alist(
    y ~ dnorm( mu , sigma ),
    mu <- a + aj[group] + (b + bj[group])*x,
    c(aj,bj)[group] ~ dmvnorm2( 0 , sigma_group , Rho_group ),
    a ~ dnorm( 0 , 10 ),
    b ~ dnorm( 0 , 1 ),
    sigma ~ dcauchy( 0 , 1 ),
    sigma_group ~ dcauchy( 0 , 1 ),
    Rho_group ~ dlkjcorr(2)
)
```


# Non-centered parameterization
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

``ensemble`` computes ``link`` and ``sim`` output for an ensemble of models, each weighted by its Akaike weight, as computed from WAIC.

