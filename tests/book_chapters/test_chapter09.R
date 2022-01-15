# test_dir('rethinking/tests/book_chapters',reporter="summary")
# test_dir('rethinking/tests/book_chapters',filter="chapter09")

context('chapter 9')
library(rethinking)

expect_equiv_eps <- function( x , y , eps=0.01 ) {
    expect_equivalent( x , y , tolerance=eps )
}


## R code 9.1
num_weeks <- 1e5
positions <- rep(0,num_weeks)
current <- 10
for ( i in 1:num_weeks ) {
    # record current position
    positions[i] <- current

    # flip coin to generate proposal
    proposal <- current + sample( c(-1,1) , size=1 )
    # now make sure he loops around the archipelago
    if ( proposal < 1 ) proposal <- 10
    if ( proposal > 10 ) proposal <- 1

    # move?
    prob_move <- proposal/current
    current <- ifelse( runif(1) < prob_move , proposal , current )
}

test_that("R code 9.1",{
    expect_equivalent( length(positions) , 1e5 )
})

## R code 9.11
library(rethinking)
data(rugged)
d <- rugged
d$log_gdp <- log(d$rgdppc_2000)
dd <- d[ complete.cases(d$rgdppc_2000) , ]
dd$log_gdp_std <- dd$log_gdp / mean(dd$log_gdp)
dd$rugged_std <- dd$rugged / max(dd$rugged)
dd$cid <- ifelse( dd$cont_africa==1 , 1 , 2 )

## R code 9.12
m8.3 <- quap(
    alist(
        log_gdp_std ~ dnorm( mu , sigma ) ,
        mu <- a[cid] + b[cid]*( rugged_std - 0.215 ) ,
        a[cid] ~ dnorm( 1 , 0.1 ) ,
        b[cid] ~ dnorm( 0 , 0.3 ) ,
        sigma ~ dexp( 1 )
    ) , data=dd )

test_that("R code 9.12",{
    expect_equivalent( dim(precis(m8.3,2)) , c(5,4) )
})

## R code 9.13
dat_slim <- list(
    log_gdp_std = dd$log_gdp_std,
    rugged_std = dd$rugged_std,
    cid = as.integer( dd$cid )
)

## R code 9.14
m9.1 <- ulam(
    alist(
        log_gdp_std ~ dnorm( mu , sigma ) ,
        mu <- a[cid] + b[cid]*( rugged_std - 0.215 ) ,
        a[cid] ~ dnorm( 1 , 0.1 ) ,
        b[cid] ~ dnorm( 0 , 0.3 ) ,
        sigma ~ dexp( 1 )
    ) , data=dat_slim , chains=1 )

test_that("R code 9.14",{
    expect_equivalent( dim(precis(m9.1,2)) , c(5,6) )
})


## R code 9.16
rm(m9.1)
m9.1 <- ulam(
    alist(
        log_gdp_std ~ dnorm( mu , sigma ) ,
        mu <- a[cid] + b[cid]*( rugged_std - 0.215 ) ,
        a[cid] ~ dnorm( 1 , 0.1 ) ,
        b[cid] ~ dnorm( 0 , 0.3 ) ,
        sigma ~ dexp( 1 )
    ) , data=dat_slim , chains=4 , cores=4 )

test_that("R code 9.16",{
    expect_equivalent( dim(precis(m9.1,2)) , c(5,6) )
})

## R code 9.22
y <- c(-1,1)
test_that("R code 9.22",{

    set.seed(11)
    m9.2 <- ulam(
        alist(
            y ~ dnorm( mu , sigma ) ,
            mu <- alpha ,
            alpha ~ dnorm( 0 , 1000 ) ,
            sigma ~ dexp( 0.0001 )
        ) , data=list(y=y) , chains=3 )

    expect_equivalent( dim(precis(m9.2,2)) , c(2,6) )
    expect_equivalent( divergent(m9.2) , 108 )
})

## R code 9.24
set.seed(11)
m9.3 <- ulam(
    alist(
        y ~ dnorm( mu , sigma ) ,
        mu <- alpha ,
        alpha ~ dnorm( 1 , 10 ) ,
        sigma ~ dexp( 1 )
    ) , data=list(y=y) , chains=3 )

test_that("R code 9.24",{
    expect_equivalent( dim(precis(m9.3,2)) , c(2,6) )
    expect_equivalent( divergent(m9.3) , 0 )
})

## R code 9.25
set.seed(41)
y <- rnorm( 100 , mean=0 , sd=1 )

## R code 9.26
test_that("R code 9.26",{
    set.seed(384)
    m9.4 <- ulam(
        alist(
            y ~ dnorm( mu , sigma ) ,
            mu <- a1 + a2 ,
            a1 ~ dnorm( 0 , 1000 ),
            a2 ~ dnorm( 0 , 1000 ),
            sigma ~ dexp( 1 )
        ) , data=list(y=y) , chains=3 )

    expect_equivalent( dim(precis(m9.4,2)) , c(3,6) )
    expect_equivalent( divergent(m9.4) , 0 )
    expect_equivalent( round(precis(m9.4,2)$n_eff[1]) , 2 )
})

## R code 9.27
test_that("R code 9.27",{

    m9.5 <- ulam(
        alist(
            y ~ dnorm( mu , sigma ) ,
            mu <- a1 + a2 ,
            a1 ~ dnorm( 0 , 10 ),
            a2 ~ dnorm( 0 , 10 ),
            sigma ~ dexp( 1 )
        ) , data=list(y=y) , chains=3 )

    expect_equivalent( dim(precis(m9.5,2)) , c(3,6) )
    expect_equivalent( divergent(m9.5) , 0 )

})

