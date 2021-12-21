# test_dir('rethinking/tests/rethinking_tests',filter="cmdstanr_threads")
# devtools::install_github("stan-dev/cmdstanr")

context('ulam_cmdstanr_threads')

library(rethinking)

N <- 1e4
x <- rnorm(N)
m <- 1 + rpois(N,2)
y <- rbinom( N , size=m , prob=inv_logit( -3 + x ) )

dat <- list( y=y , x=x , m=m )

# single thread
m1 <- ulam(
    alist(
        y ~ binomial_logit( m , logit_p ),
        logit_p <- a + b*x,
        a ~ normal(0,1.5),
        b ~ normal(0,0.5)
    ) , data=dat , 
    cmdstan=TRUE , threads=1 , chains=1 , cores=1 , refresh=1000 )

# two threads
m2 <- ulam(
    alist(
        y ~ binomial_logit( m , logit_p ),
        logit_p <- a + b*x,
        a ~ normal(0,1.5),
        b ~ normal(0,0.5)
    ) , data=dat , 
    cmdstan=TRUE , threads=2 , chains=2 , cores=2 , refresh=1000 )

# hierarchical version
N_id <- round(N/10)
id <- rep(1:N_id,each=10)
a <- rnorm(N_id,0,1.5)
y <- rbinom( N , size=m , prob=inv_logit( -3 + a[id] + x ) )
dat$id <- id
dat$y <- y

system.time(
m3 <- ulam(
    alist(
        y ~ binomial_logit( m , logit_p ),
        logit_p <- a + z[id]*tau + b*x,
        a ~ normal(0,1.5),
        b ~ normal(0,0.5),
        z[id] ~ normal(0,1),
        tau ~ exponential(1)
    ) , data=dat , 
    cmdstan=TRUE , threads=2 , chains=1 , cores=1 , refresh=1000 ) )

###########
# ordered logistic
# problems with initialization (just cmdstan not threading)
# problems on gcc

library(rethinking)
data(Trolley)
d <- Trolley

dat <- list(
    R = d$response,
    A = d$action,
    I = d$intention,
    C = d$contact )

system.time(
m12.5 <- ulam(
    alist(
        R ~ dordlogit( phi , cutpoints ),
        phi <- bA*A + bC*C + (bI + bIA*A + bIC*C)*I ,
        c(bA,bI,bC,bIA,bIC) ~ dnorm( 0 , 0.5 ),
        cutpoints ~ dnorm( 0 , 1.5 )
    ) , data=dat , chains=1 , cores=1 , cmdstan=TRUE , threads=4 , sample=TRUE , refresh=500 , grain=1 ) )

m0 <- ulam(
    alist(
        y ~ dordlogit(0,cutpoints),
        ordered[6]:cutpoints ~ normal(0,1.5)
    ) , data=list(y=1:4), cmdstan=TRUE , sample=TRUE )

########
# poisson GLM - m11.11

# fixed link() method N issue

library(rethinking)
data(Kline)
d <- Kline
d$contact_id <- ifelse( d$contact=="high" , 2 , 1 )
dat2 <- list( T=d$total_tools, P=as.numeric(d$population), cid=d$contact_id )

m0 <- ulam(
    alist(
        T ~ dpois( lambda ),
        lambda <- exp(a[cid])*P^b[cid]/g,
        a[cid] ~ dnorm(1,1),
        b[cid] ~ dexp(1),
        g ~ dexp(1)
    ), data=dat2 , chains=4 , cores=4 , cmdstan=TRUE , threads=2 , sample=TRUE , refresh=1000 )

precis(m0,3)
l <- link(m0)

## R code 14.37
# load the distance matrix
library(rethinking)
data(islandsDistMatrix)

# display (measured in thousands of km)
Dmat <- islandsDistMatrix
colnames(Dmat) <- c("Ml","Ti","SC","Ya","Fi","Tr","Ch","Mn","To","Ha")

## R code 14.39
data(Kline2) # load the ordinary data, now with coordinates
d <- Kline2
d$society <- 1:10 # index observations

dat_list <- list(
    T = d$total_tools,
    P = d$population,
    society = d$society,
    Dmat=islandsDistMatrix )

# set_ulam_cmdstan( TRUE )

# have to assign SIGMA to transpars to get this to work
m1 <- ulam(
    alist(
        T ~ dpois(lambda),
        lambda <- (a*P^b/g)*exp(k[society]),
        vector[10]:k ~ multi_normal( 0 , SIGMA ),
        transpars> matrix[10,10]:SIGMA <- cov_GPL2( Dmat , etasq , rhosq , 0.01 ),
        c(a,b,g) ~ dexp( 1 ),
        etasq ~ dexp( 2 ),
        rhosq ~ dexp( 0.5 )
    ), data=dat_list , chains=1 , cores=1 , iter=2000 , threads=2 , sample=TRUE , cmdstan=TRUE )


######
# varying intercepts

library(rethinking)
data(chimpanzees)
d <- chimpanzees
d$treatment <- 1 + d$prosoc_left + 2*d$condition

dat_list <- list(
    pulled_left = d$pulled_left,
    actor = d$actor,
    block_id = d$block,
    treatment = as.integer(d$treatment) )

m13.4 <- ulam(
    alist(
        pulled_left ~ dbinom( 1 , p ) ,
        logit(p) <- a[actor] + g[block_id] + b[treatment] ,
        b[treatment] ~ dnorm( 0 , 0.5 ),
      ## adaptive priors
        a[actor] ~ dnorm( a_bar , sigma_a ),
        g[block_id] ~ dnorm( 0 , sigma_g ),
      ## hyper-priors
        a_bar ~ dnorm( 0 , 1.5 ),
        sigma_a ~ dexp(1),
        sigma_g ~ dexp(1)
    ) , data=dat_list , chains=1 , cores=1 , cmdstan=TRUE , threads=2 )

library(rethinking)
data(chimpanzees)
d <- chimpanzees
d$block_id <- d$block
d$treatment <- 1L + d$prosoc_left + 2L*d$condition
dat <- list(
    L = d$pulled_left,
    tid = d$treatment,
    actor = d$actor,
    block_id = as.integer(d$block_id) )

m14.3 <- ulam(
    alist(
        L ~ binomial(1,p),
        logit(p) <- g[tid] + alpha[actor,tid] + beta[block_id,tid],

        # adaptive priors - non-centered
        transpars> matrix[actor,4]:alpha <-
                compose_noncentered( sigma_actor , L_Rho_actor , z_actor ),
        transpars> matrix[block_id,4]:beta <-
                compose_noncentered( sigma_block , L_Rho_block , z_block ),
        matrix[4,actor]:z_actor ~ normal( 0 , 1 ),
        matrix[4,block_id]:z_block ~ normal( 0 , 1 ),

        # fixed priors
        g[tid] ~ normal(0,1),
        vector[4]:sigma_actor ~ dexp(1),
        cholesky_factor_corr[4]:L_Rho_actor ~ lkj_corr_cholesky( 2 ),
        vector[4]:sigma_block ~ dexp(1),
        cholesky_factor_corr[4]:L_Rho_block ~ lkj_corr_cholesky( 2 ),

        # compute ordinary correlation matrixes from Cholesky factors
        gq> matrix[4,4]:Rho_actor <<- Chol_to_Corr(L_Rho_actor),
        gq> matrix[4,4]:Rho_block <<- Chol_to_Corr(L_Rho_block)
    ) , data=dat , chains=4 , cores=4 , cmdstan=TRUE , threads=2 , sample=TRUE )


