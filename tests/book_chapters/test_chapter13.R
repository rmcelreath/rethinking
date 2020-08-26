# test_dir('rethinking/tests/book_chapters',reporter="summary")
# test_dir('rethinking/tests/book_chapters',filter="chapter13")

context('chapter 13')
library(rethinking)

expect_equiv_eps <- function( x , y , eps=0.01 ) {
    expect_equivalent( x , y , tolerance=eps )
}

# set_ulam_cmdstan(TRUE)

## R code 13.1
library(rethinking)
data(reedfrogs)
d <- reedfrogs

## R code 13.2
# make the tank cluster variable
d$tank <- 1:nrow(d)

dat <- list(
    S = d$surv,
    N = d$density,
    tank = d$tank )

# approximate posterior
m13.1 <- ulam(
    alist(
        S ~ dbinom( N , p ) ,
        logit(p) <- a[tank] ,
        a[tank] ~ dnorm( 0 , 1.5 )
    ), data=dat , chains=4 , log_lik=TRUE )

test_that("R code 13.2",{
    expect_equivalent( dim(precis(m13.1,2)) , c(48,6) )
})

## R code 13.3
m13.2 <- ulam(
    alist(
        S ~ dbinom( N , p ) ,
        logit(p) <- a[tank] ,
        a[tank] ~ dnorm( a_bar , sigma ) ,
        a_bar ~ dnorm( 0 , 1.5 ) ,
        sigma ~ dexp( 1 )
    ), data=dat , chains=4 , log_lik=TRUE )

test_that("R code 13.3",{
    expect_equivalent( dim(precis(m13.2,2)) , c(50,6) )
})

## R code 13.4
test_that("R code 13.4",{
    xxx <- compare( m13.1 , m13.2 )
    expect_equivalent( dim(xxx) , c(2,6) )
})

## R code 13.21
library(rethinking)
data(chimpanzees)
d <- chimpanzees
d$treatment <- 1 + d$prosoc_left + 2*d$condition

dat_list <- list(
    pulled_left = d$pulled_left,
    actor = d$actor,
    block_id = d$block,
    treatment = as.integer(d$treatment) )

test_that("R code 13.21",{

    set.seed(13)
    expect_warning( 
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
        ) , data=dat_list , chains=4 , cores=4 , log_lik=TRUE )
    )

    expect_equivalent( dim(precis(m13.4,2)) , c(20,6) )
})

## R code 13.23

test_that("R code 13.23",{

    set.seed(14)
    m13.5 <- ulam(
        alist(
            pulled_left ~ dbinom( 1 , p ) ,
            logit(p) <- a[actor] + b[treatment] ,
            b[treatment] ~ dnorm( 0 , 0.5 ),
            a[actor] ~ dnorm( a_bar , sigma_a ),
            a_bar ~ dnorm( 0 , 1.5 ),
            sigma_a ~ dexp(1)
        ) , data=dat_list , chains=4 , cores=4 , log_lik=TRUE )

    expect_equivalent( dim(precis(m13.5,2)) , c(13,6) )

})

## R code 13.25

test_that("R code 13.25",{

    set.seed(15)
    expect_warning(
        m13.6 <- ulam(
            alist(
                pulled_left ~ dbinom( 1 , p ) ,
                logit(p) <- a[actor] + g[block_id] + b[treatment] ,
                b[treatment] ~ dnorm( 0 , sigma_b ),
                a[actor] ~ dnorm( a_bar , sigma_a ),
                g[block_id] ~ dnorm( 0 , sigma_g ),
                a_bar ~ dnorm( 0 , 1.5 ),
                sigma_a ~ dexp(1),
                sigma_g ~ dexp(1),
                sigma_b ~ dexp(1)
            ) , data=dat_list , chains=4 , cores=4 , log_lik=TRUE )
    )

    expect_equivalent( dim(precis(m13.6,2)) , c(21,6) )
})

## R code 13.26

test_that("R code 13.26",{
    expect_warning(
        m13.7 <- ulam(
            alist(
                v ~ normal(0,3),
                x ~ normal(0,exp(v))
            ), data=list(N=1) , chains=4 )
    )
    expect_equivalent( dim(precis(m13.7,2)) , c(2,6) )
})

## R code 13.27

test_that("R code 13.27",{
    m13.7nc <- ulam(
        alist(
            v ~ normal(0,3),
            z ~ normal(0,1),
            gq> real[1]:x <<- z*exp(v)
        ), data=list(N=1) , chains=4 )
    expect_equivalent( dim(precis(m13.7nc,2)) , c(3,6) )
})

## R code 13.28

test_that("R code 13.28",{
    set.seed(13)
    expect_warning(
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
            ) , data=dat_list , chains=4 , cores=4 , log_lik=TRUE )
    )
    expect_warning(
        m13.4b <- ulam( m13.4 , chains=4 , cores=4 , control=list(adapt_delta=0.99) )
    )

    expect_equivalent( dim(precis(m13.4b,2)) , c(20,6) )
})

## R code 13.29

test_that("R code 13.29",{

    set.seed(13)
    #expect_warning(
        m13.4nc <- ulam(
            alist(
                pulled_left ~ dbinom( 1 , p ) ,
                logit(p) <- a_bar + z[actor]*sigma_a + # actor intercepts
                            x[block_id]*sigma_g +      # block intercepts
                            b[treatment] ,
                b[treatment] ~ dnorm( 0 , 0.5 ),
                z[actor] ~ dnorm( 0 , 1 ),
                x[block_id] ~ dnorm( 0 , 1 ),
                a_bar ~ dnorm( 0 , 1.5 ),
                sigma_a ~ dexp(1),
                sigma_g ~ dexp(1),
                gq> vector[actor]:a <<- a_bar + z*sigma_a,
                gq> vector[block_id]:g <<- x*sigma_g
            ) , data=dat_list , chains=4 , cores=4 )
    #)

    expect_equivalent( dim(precis(m13.4nc,2)) , c(33,6) )

})


