# test --O1


N <- 1e4
x <- rnorm(N)
m <- 1 + rpois(N,2)
y <- rbinom( N , size=m , prob=inv_logit( -3 + x ) )
dat <- list( y=y , x=x , m=m )

N_id <- round(N/10)
id <- rep(1:N_id,each=10)
a <- rnorm(N_id,0,1.5)
y <- rbinom( N , size=m , prob=inv_logit( -3 + a[id] + x ) )
dat$id <- id
dat$y <- y

system.time(
mO1 <- ulam(
    alist(
        y ~ binomial_logit( m , logit_p ),
        logit_p <- a + z[id]*tau + b*x,
        a ~ normal(0,1.5),
        b ~ normal(0,0.5),
        z[id] ~ normal(0,1),
        tau ~ exponential(1)
    ) , data=dat , 
    cmdstan=TRUE , chains=4 , threads=1 , cores=4 , refresh=1000 , stanc_options = list("O1") ) )

system.time(
mO0 <- ulam(
    alist(
        y ~ binomial_logit( m , logit_p ),
        logit_p <- a + z[id]*tau + b*x,
        a ~ normal(0,1.5),
        b ~ normal(0,0.5),
        z[id] ~ normal(0,1),
        tau ~ exponential(1)
    ) , data=dat , 
    cmdstan=TRUE , chains=4 , threads=1 , cores=4 , refresh=1000 , stanc_options = list("O0") ) )

# slopes

dat$N_id <- N_id
system.time(
mO1s <- ulam(
    alist(
        y ~ binomial_logit( m , logit_p ),
        logit_p <- a + v[id,1] + (b+v[id,2])*x,
        transpars> matrix[N_id,2]:v <- compose_noncentered( tau , LR , z ),
        a ~ normal(0,1.5),
        b ~ normal(0,0.5),
        matrix[2,N_id]:z ~ normal(0,1),
        vector[2]:tau ~ dexp(1),
        cholesky_factor_corr[2]:LR ~ lkj_corr_cholesky( 2 )
    ) , data=dat , 
    cmdstan=TRUE , chains=4 , threads=1 , cores=4 , refresh=1000 , stanc_options = list("O1") ) )

system.time(
mO0s <- ulam(
    alist(
        y ~ binomial_logit( m , logit_p ),
        logit_p <- a + v[id,1] + (b+v[id,2])*x,
        transpars> matrix[N_id,2]:v <- compose_noncentered( tau , LR , z ),
        a ~ normal(0,1.5),
        b ~ normal(0,0.5),
        matrix[2,N_id]:z ~ normal(0,1),
        vector[2]:tau ~ dexp(1),
        cholesky_factor_corr[2]:LR ~ lkj_corr_cholesky( 2 )
    ) , data=dat , 
    cmdstan=TRUE , chains=4 , threads=1 , cores=4 , refresh=1000 , stanc_options = list("O0") ) )

