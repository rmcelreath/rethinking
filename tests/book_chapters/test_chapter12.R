# test_dir('rethinking/tests/book_chapters',reporter="summary")
# test_dir('rethinking/tests/book_chapters',filter="chapter12")

context('chapter 12')
library(rethinking)

expect_equiv_eps <- function( x , y , eps=0.01 ) {
    expect_equivalent( x , y , tolerance=eps )
}


## R code 12.2
library(rethinking)
data(UCBadmit)
d <- UCBadmit
d$gid <- ifelse( d$applicant.gender=="male" , 1L , 2L )
dat <- list( A=d$admit , N=d$applications , gid=d$gid )
m12.1 <- ulam(
    alist(
        A ~ dbetabinom( N , pbar , theta ),
        logit(pbar) <- a[gid],
        a[gid] ~ dnorm( 0 , 1.5 ),
        transpars> theta <<- phi + 2.0,
        phi ~ dexp(1)
    ), data=dat , chains=4 )

test_that("R code 12.2",{
    expect_equivalent( dim(precis(m12.1,2)) , c(4,6) )
})

## R code 12.5
# postcheck( m12.1 )

## R code 12.6
library(rethinking)
data(Kline)
d <- Kline
d$P <- standardize( log(d$population) )
d$contact_id <- ifelse( d$contact=="high" , 2L , 1L )

dat2 <- list(
    T = d$total_tools,
    P = d$population,
    cid = d$contact_id )

m12.2 <- ulam(
    alist(
        T ~ dgampois( lambda , phi ),
        lambda <- exp(a[cid])*P^b[cid] / g,
        a[cid] ~ dnorm(1,1),
        b[cid] ~ dexp(1),
        g ~ dexp(1),
        phi ~ dexp(1)
    ), data=dat2 , chains=4 , log_lik=TRUE )

test_that("R code 12.6",{
    expect_equivalent( dim(precis(m12.2,2)) , c(6,6) )
})

## R code 12.7
# define parameters
prob_drink <- 0.2 # 20% of days
rate_work <- 1    # average 1 manuscript per day

# sample one year of production
N <- 365

# simulate days monks drink
set.seed(365)
drink <- rbinom( N , 1 , prob_drink )

# simulate manuscripts completed
y <- (1-drink)*rpois( N , rate_work )

## R code 12.9
m12.3 <- ulam(
    alist(
        y ~ dzipois( p , lambda ),
        logit(p) <- ap,
        log(lambda) <- al,
        ap ~ dnorm( -1.5 , 1 ),
        al ~ dnorm( 1 , 0.5 )
    ) , data=list(y=y) , chains=4 )

test_that("R code 12.9",{
    expect_equivalent( dim(precis(m12.3,2)) , c(2,6) )
})

## R code 12.11
m12.3_alt <- ulam(
    alist(
        y|y>0 ~ custom( log1m(p) + poisson_lpmf(y|lambda) ),
        y|y==0 ~ custom( log_mix( p , 0 , poisson_lpmf(0|lambda) ) ),
        logit(p) <- ap,
        log(lambda) <- al,
        ap ~ dnorm(-1.5,1),
        al ~ dnorm(1,0.5)
    ) , data=list(y=as.integer(y)) , chains=4 )

test_that("R code 12.9",{
    expect_equivalent( dim(precis(m12.3_alt,2)) , c(2,6) )
})

## R code 12.12
library(rethinking)
data(Trolley)
d <- Trolley

## R code 12.16
m12.4 <- ulam(
    alist(
        R ~ dordlogit( 0 , cutpoints ),
        cutpoints ~ dnorm( 0 , 1.5 )
    ) , data=list( R=d$response ), chains=4 , cores=4 )

test_that("R code 12.16",{
    expect_equivalent( dim(precis(m12.4,2)) , c(6,6) )
})

## R code 12.17
m12.4q <- quap(
    alist(
        response ~ dordlogit( 0 , c(a1,a2,a3,a4,a5,a6) ),
        c(a1,a2,a3,a4,a5,a6) ~ dnorm( 0 , 1.5 )
    ) , data=d , start=list(a1=-2,a2=-1,a3=0,a4=1,a5=2,a6=2.5) )

test_that("R code 12.17",{
    expect_equivalent( dim(precis(m12.4q,2)) , c(6,4) )
})

## R code 12.24
dat <- list(
    R = d$response,
    A = d$action,
    I = d$intention,
    C = d$contact )
m12.5 <- ulam(
    alist(
        R ~ dordlogit( phi , cutpoints ),
        phi <- bA*A + bC*C + BI*I ,
        BI <- bI + bIA*A + bIC*C ,
        c(bA,bI,bC,bIA,bIC) ~ dnorm( 0 , 0.5 ),
        cutpoints ~ dnorm( 0 , 1.5 )
    ) , data=dat , chains=4 , cores=4 )

test_that("R code 12.24",{
    expect_equivalent( dim(precis(m12.5,2)) , c(11,6) )
})

## R code 12.30
library(rethinking)
data(Trolley)
d <- Trolley

## R code 12.31
edu_levels <- c( 6 , 1 , 8 , 4 , 7 , 2 , 5 , 3 )
d$edu_new <- edu_levels[ d$edu ]

## R code 12.34
dat <- list(
    R = d$response ,
    action = d$action,
    intention = d$intention,
    contact = d$contact,
    E = as.integer( d$edu_new ),   # edu_new as an index
    alpha = rep( 2 , 7 ) )      # delta prior

m12.6 <- ulam(
    alist(
        R ~ ordered_logistic( phi , kappa ),
        phi <- bE*sum( delta_j[1:E] ) + bA*action + bI*intention + bC*contact,
        kappa ~ normal( 0 , 1.5 ),
        c(bA,bI,bC,bE) ~ normal( 0 , 1 ),
        vector[8]: delta_j <<- append_row( 0 , delta ),
        simplex[7]: delta ~ dirichlet( alpha )
    ), data=dat , chains=4 , cores=4 )

test_that("R code 12.34",{
    expect_equivalent( dim(precis(m12.6,2)) , c(17,6) )
})

## R code 12.37
dat$edu_norm <- normalize( d$edu_new )
m12.7 <- ulam(
    alist(
        R ~ ordered_logistic( mu , cutpoints ),
        mu <- bE*edu_norm + bA*action + bI*intention + bC*contact,
        c(bA,bI,bC,bE) ~ normal( 0 , 1 ),
        cutpoints ~ normal( 0 , 1.5 )
    ), data=dat , chains=4 , cores=4 )

test_that("R code 12.37",{
    expect_equivalent( dim(precis(m12.7,2)) , c(10,6) )
})


