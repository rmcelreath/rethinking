wythamewa
==========

This package contains models and data to replicate the payoff and frequency dependent social learning analysis of the Aplin, Sheldon, McElreath 2017, "Conformity does not perpetuate suboptimal traditions in a wild population of songbirds" (doi: 10.1073/pnas.1621067114).

This package requires ``rstan`` (http://mc-stan.org/) and ``rethinking`` (https://github.com/rmcelreath/rethinking) to fit the statistical models.

See ``?wythamewa`` for a quick start guide to replicating the analysis and figures in the paper. Full scripts to reproduce the figures are found in the various help pages accessible from there.

The full model can be replicated with the following code:
```
lms <- list(
    "mu[1] + a_bird[bird[i],1] + b_age[1]*age[i]",
    "mu[2] + a_bird[bird[i],2] + b_age[2]*age[i]",
    "mu[3] + a_bird[bird[i],3] + b_age[3]*age[i]",
    "mu[4] + a_bird[bird[i],4] + b_age[4]*age[i]"
    )
links <- c("logit", "logit", "log", "logit", "")
prior <- "
    mu ~ normal(0,1);
    diff_hi ~ cauchy(0,1);
    b_age ~ normal(0,1);
    to_vector(z_bird) ~ normal(0,1);
    L_Rho_bird ~ lkj_corr_cholesky(3);
    sigma_bird ~ exponential(2);"
mod1 <- ewa_def( model=lms , prior=prior , link=links )
set.seed(1)
m <- ewa_fit( mod1 , warmup=500 , iter=1000 , chains=3 , cores=3 ,
    control=list( adapt_delta=0.99 , max_treedepth=12 ) )
```
