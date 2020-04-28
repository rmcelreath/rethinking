# test_dir('rethinking/tests/testthat')
# test_dir('rethinking/tests/book_chapters',filter="chapter04", reporter="summary")

context('chapter 4')
library(rethinking)

expect_equiv_eps <- function( x , y , eps=0.01 ) {
    expect_equivalent( x , y , tolerance=eps )
}

## R code 4.26
library(rethinking)
data(Howell1)
d <- Howell1
d2 <- d[ d$age >= 18 , ]

## R code 4.27
flist <- alist(
    height ~ dnorm( mu , sigma ) ,
    mu ~ dnorm( 178 , 20 ) ,
    sigma ~ dunif( 0 , 50 )
)

## R code 4.28
m4.1 <- quap( flist , data=d2 )

## R code 4.29
test_that("R code 4.29",
    expect_equiv_eps( round(coef(m4.1),2) , c(154.61,7.73) )
)

## R code 4.31
m4.2 <- quap(
    alist(
        height ~ dnorm( mu , sigma ) ,
        mu ~ dnorm( 178 , 0.1 ) ,
        sigma ~ dunif( 0 , 50 )
    ) , data=d2 )

test_that("R code 4.31",
    expect_equiv_eps( round(coef(m4.2),2) , c(177.86,24.52) )
)

## R code 4.42
# load data again, since it's a long way back
library(rethinking)
data(Howell1)
d <- Howell1
d2 <- d[ d$age >= 18 , ]

# define the average weight, x-bar
d2$xbar <- mean(d2$weight)

# fit model
m4.3 <- quap(
    alist(
        height ~ dnorm( mu , sigma ) ,
        mu <- a + b*( weight - xbar ) ,
        a ~ dnorm( 178 , 20 ) ,
        b ~ dlnorm( 0 , 1 ) ,
        sigma ~ dunif( 0 , 50 )
    ) ,
    data=d2 )

test_that("R code 4.42",
    expect_equiv_eps( round(coef(m4.3),2) , c(154.60  , 0.90  , 5.07) )
)

## R code 4.43
m4.3b <- quap(
    alist(
        height ~ dnorm( mu , sigma ) ,
        mu <- a + exp(log_b)*( weight - xbar ),
        a ~ dnorm( 178 , 20 ) ,
        log_b ~ dnorm( 0 , 1 ) ,
        sigma ~ dunif( 0 , 50 )
    ) ,
    data=d2 )

test_that("R code 4.43",
    expect_equiv_eps( round(coef(m4.3b),2) , c(154.60  , -0.10  , 5.07) )
)

## R code 4.64
library(rethinking)
data(Howell1)
d <- Howell1

## R code 4.65
d$weight_s <- ( d$weight - mean(d$weight) )/sd(d$weight)
d$weight_s2 <- d$weight_s^2
m4.5 <- quap(
    alist(
        height ~ dnorm( mu , sigma ) ,
        mu <- a + b1*weight_s + b2*weight_s2 ,
        a ~ dnorm( 178 , 20 ) ,
        b1 ~ dlnorm( 0 , 1 ) ,
        b2 ~ dnorm( 0 , 1 ) ,
        sigma ~ dunif( 0 , 50 )
    ) ,
    data=d )

test_that("R code 4.65",
    expect_equiv_eps( round(coef(m4.5),2) , 
    c(146.06 , 21.73 , -7.80  , 5.77) )
)

## R code 4.69
d$weight_s3 <- d$weight_s^3
m4.6 <- quap(
    alist(
        height ~ dnorm( mu , sigma ) ,
        mu <- a + b1*weight_s + b2*weight_s2 + b3*weight_s3 ,
        a ~ dnorm( 178 , 20 ) ,
        b1 ~ dlnorm( 0 , 1 ) ,
        b2 ~ dnorm( 0 , 10 ) ,
        b3 ~ dnorm( 0 , 10 ) ,
        sigma ~ dunif( 0 , 50 )
    ) ,
    data=d )

test_that("R code 4.69",
    expect_equiv_eps( round(coef(m4.6),2) , 
    c(146.74 , 15.00 , -6.54 ,  3.60 ,  4.82) )
)

## R code 4.72
library(rethinking)
data(cherry_blossoms)
d <- cherry_blossoms

## R code 4.73
d2 <- d[ complete.cases(d$temp) , ] # complete cases on temp
num_knots <- 15
knot_list <- quantile( d2$year , probs=seq(0,1,length.out=num_knots) )

## R code 4.74
library(splines)
B <- bs(d2$year,
    knots=knot_list[-c(1,num_knots)] ,
    degree=3 , intercept=TRUE )

## R code 4.76
m4.7 <- quap(
    alist(
        T ~ dnorm( mu , sigma ) ,
        mu <- a + B %*% w ,
        a ~ dnorm(6,10),
        w ~ dnorm(0,1),
        sigma ~ dexp(1)
    ),
    data=list( T=d2$temp , B=B ) ,
    start=list( w=rep( 0 , ncol(B) ) ) )

test_that("R code 4.76",
    expect_equiv_eps( round(coef(m4.7),2) , 
    c(0.10 , 0.24 , 1.20, -0.82 , 0.10, -1.40 , 1.16, -1.92 , 2.35 ,-2.32 , 0.94, -1.61 , 0.18 , -1.24 , 0.03 , 1.00 , 2.04 , 6.32 , 0.34 ) 
    )
)

## R code 4.79
m4.7alt <- quap(
    alist(
        T ~ dnorm( mu , sigma ) ,
        mu <- a + sapply( 1:1124 , function(i) sum( B[i,]*w ) ) ,
        a ~ dnorm(6,10),
        w ~ dnorm(0,1),
        sigma ~ dexp(1)
    ),
    data=list( T=d2$temp , B=B ) ,
    start=list( w=rep( 0 , ncol(B) ) ) )

test_that("R code 4.76",
    expect_equiv_eps( round(coef(m4.7alt),2) , 
    c(0.10 , 0.24 , 1.20, -0.82 , 0.10, -1.40 , 1.16, -1.92 , 2.35 ,-2.32 , 0.94, -1.61 , 0.18 , -1.24 , 0.03 , 1.00 , 2.04 , 6.32 , 0.34 ) 
    )
)
