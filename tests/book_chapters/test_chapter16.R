# test_dir('rethinking/tests/book_chapters',reporter="summary")
# test_dir('rethinking/tests/book_chapters',filter="chapter16")

context('chapter 16')
library(rethinking)

expect_equiv_eps <- function( x , y , eps=0.01 ) {
    expect_equivalent( x , y , tolerance=eps )
}


## R code 16.1
library(rethinking)
data(Howell1)
d <- Howell1

# scale observed variables
d$w <- d$weight / mean(d$weight)
d$h <- d$height / mean(d$height)

## R code 16.2
test_that("R code 16.2",{
    m16.1 <- ulam(
        alist(
            w ~ dlnorm( mu , sigma ),
            exp(mu) <- 3.141593 * k * p^2 * h^3,
            p ~ beta( 2 , 18 ),
            k ~ exponential( 0.5 ),
            sigma ~ exponential( 1 )
        ), data=d , chains=4 , cores=4 )
    expect_equivalent( dim(precis(m16.1,3)) , c(3,6) )
})

## R code 16.4
library(rethinking)
data(Boxes)

## R code 16.7
data(Boxes_model)

## R code 16.8
# prep data
dat_list <- list(
    N = nrow(Boxes),
    y = Boxes$y,
    majority_first = Boxes$majority_first )

# run the sampler
test_that("R code 16.8",{
    # m16.2 <- stan( model_code=Boxes_model , data=dat_list , chains=3 , cores=3 )
    m16.2 <- cstan( model_code=Boxes_model , data=dat_list , chains=3 , cores=3 )
    expect_equivalent( dim(precis(m16.2,3)) , c(5,6) )
})

## R code 16.9
library(rethinking)
data(Panda_nuts)

## R code 16.11
dat_list <- list(
    n = as.integer( Panda_nuts$nuts_opened ),
    age = Panda_nuts$age / max(Panda_nuts$age),
    seconds = Panda_nuts$seconds )

test_that("R code 16.11",{
    m16.4 <- ulam(
        alist(
            n ~ poisson( lambda ),
            lambda <- seconds*phi*(1-exp(-k*age))^theta,
            phi ~ lognormal( log(1) , 0.1 ),
            k ~ lognormal( log(2) , 0.25 ),
            theta ~ lognormal( log(5) , 0.25 )
        ), data=dat_list , chains=4 )
    expect_equivalent( dim(precis(m16.4,3)) , c(3,6) )
})

## R code 16.13
library(rethinking)
data(Lynx_Hare)

## R code 16.17
data(Lynx_Hare_model)

## R code 16.18
dat_list <- list(
    N = nrow(Lynx_Hare),
    pelts = Lynx_Hare[,2:3] )

test_that("R code 16.18",{
    m16.5 <- cstan( model_code=Lynx_Hare_model , data=dat_list , chains=3 , cores=3 , control=list( adapt_delta=0.95 ) )
    expect_equivalent( dim(precis(m16.5,3)) , c(94,6) )
})

