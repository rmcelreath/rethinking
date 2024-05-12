# test_dir('rethinking/tests/rethinking_tests',filter="cmdstanr_threads")
# devtools::install_github("stan-dev/cmdstanr")

context('cstan_cmdstanr_support')

library(rethinking)

N <- 1e4
x <- rnorm(N)
m <- 1 + rpois(N,2)
y <- rbinom( N , size=m , prob=inv_logit( -3 + x ) )

dat <- list( y=y , x=x , m=m )

# model code
mc <- "
data{
    array[10000] int m;
    array[10000] int y;
     vector[10000] x;
}
parameters{
     real a;
     real b;
     vector[4] theta;
}
model{
     vector[10000] logit_p;
    theta ~ normal(0,1);
    b ~ normal( 0 , 0.5 );
    a ~ normal( 0 , 1.5 );
    for ( i in 1:10000 ) {
        logit_p[i] = a + b * x[i];
    }
    y ~ binomial_logit( m , logit_p );
}
"

m <- cstan( model_code=mc , data=dat , chains=4 , cores=4 , rstan_out=FALSE )

precis(m,2)

post <- extract.samples(m)

# with start inits

m2 <- cstan( model_code=mc , data=dat , chains=4 , cores=4 , rstan_out=FALSE , start=list(a=0,b=0,theta=rep(0,4)) )

precis(m2,2)


