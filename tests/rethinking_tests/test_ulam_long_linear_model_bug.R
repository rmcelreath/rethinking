# test for ulam long linear model bug

library(rethinking)

N <- 100
x <- rnorm(N)
y <- rnorm(N,x)

m <- ulam(
    alist(
        y ~ normal(mu,sigma),
        mu <- 
            intercept + beta_slope_x*x + 
            intercept + beta_slope_x*x + 
            intercept + beta_slope_x*x + 
            intercept + beta_slope_x*x + 
            intercept + beta_slope_x*x + 
            intercept + beta_slope_x*x + 
            intercept + beta_slope_x*x + 
            intercept + beta_slope_x*x + 
            intercept + beta_slope_x*x + 
            intercept + beta_slope_x*x + 
            intercept + beta_slope_x*x + 
            intercept + beta_slope_x*x + 
            intercept + beta_slope_x*x + 
            intercept + beta_slope_x*x + 
            intercept + beta_slope_x*x + 
            intercept + beta_slope_x*x + 
            intercept + beta_slope_x*x + 
            intercept + beta_slope_x*x + 
            intercept + beta_slope_x*x + 
            intercept + beta_slope_x*x,
        intercept ~ normal(0,0.5),
        beta_slope_x ~ normal(0,0.5),
        sigma ~ exponential(1)
    ), data=list(x=x,y=y) , chains=1 )
