# test_dir('rethinking/tests/book_chapters',reporter="summary")
# test_dir('rethinking/tests/book_chapters',filter="chapter05")

context('chapter 5')
library(rethinking)
library(openssl) # for md5 hash

expect_equiv_eps <- function( x , y , eps=0.01 ) {
    expect_equivalent( x , y , tolerance=eps )
}

## R code 5.1
# load data and copy
library(rethinking)
data(WaffleDivorce)
d <- WaffleDivorce

# standardize variables
d$A <- scale( d$MedianAgeMarriage )
d$D <- scale( d$Divorce )

## R code 5.3
m5.1 <- quap(
    alist(
        D ~ dnorm( mu , sigma ) ,
        mu <- a + bA * A ,
        a ~ dnorm( 0 , 0.2 ) ,
        bA ~ dnorm( 0 , 0.5 ) ,
        sigma ~ dexp( 1 )
    ) , data = d )

test_that("R code 5.3",
    expect_equiv_eps( round(coef(m5.1),2) , c( 0.00 ,-0.57 , 0.79 ) )
)

## R code 5.6
d$M <- scale( d$Marriage )
m5.2 <- quap(
    alist(
        D ~ dnorm( mu , sigma ) ,
        mu <- a + bM * M ,
        a ~ dnorm( 0 , 0.2 ) ,
        bM ~ dnorm( 0 , 0.5 ) ,
        sigma ~ dexp( 1 )
    ) , data = d )

test_that("R code 5.6",
    expect_equiv_eps( round(coef(m5.2),2) , c(  0.00 , 0.35 , 0.91 ) )
)

## R code 5.10
m5.3 <- quap(
    alist(
        D ~ dnorm( mu , sigma ) ,
        mu <- a + bM*M + bA*A ,
        a ~ dnorm( 0 , 0.2 ) ,
        bM ~ dnorm( 0 , 0.5 ) ,
        bA ~ dnorm( 0 , 0.5 ) ,
        sigma ~ dexp( 1 )
    ) , data = d )

test_that("R code 5.10",
    expect_equiv_eps( round(coef(m5.3),2) , c( 0.00, -0.07, -0.61 , 0.79) )
)

## R code 5.13
m5.4 <- quap(
    alist(
        M ~ dnorm( mu , sigma ) ,
        mu <- a + bAM * A ,
        a ~ dnorm( 0 , 0.2 ) ,
        bAM ~ dnorm( 0 , 0.5 ) ,
        sigma ~ dexp( 1 )
    ) , data = d )

## R code 5.14
mu <- link(m5.4)
mu_mean <- apply( mu , 2 , mean )
mu_resid <- d$M - mu_mean

## R code 5.15
# call link without specifying new data
# so it uses original data
mu <- link( m5.3 )

# summarize samples across cases
mu_mean <- apply( mu , 2 , mean )
mu_PI <- apply( mu , 2 , PI )

# simulate observations
# again no new data, so uses original data
D_sim <- sim( m5.3 , n=1e4 )
D_PI <- apply( D_sim , 2 , PI )

test_that("R code 5.15", {
    expect_equivalent( dim(mu) , c(1000,50) )
    expect_equivalent( length(mu_mean) , 50 )
    expect_equivalent( dim(mu_PI) , c(2,50) )
    expect_equivalent( dim(D_sim) , c(10000,50) )
    expect_equivalent( dim(D_PI) , c(2,50) )
})

## R code 5.19
data(WaffleDivorce)
d <- list()
d$A <- standardize( WaffleDivorce$MedianAgeMarriage )
d$D <- standardize( WaffleDivorce$Divorce )
d$M <- standardize( WaffleDivorce$Marriage )

m5.3_A <- quap(
    alist(
      ## A -> D <- M
        D ~ dnorm( mu , sigma ) ,
        mu <- a + bM*M + bA*A ,
        a ~ dnorm( 0 , 0.2 ) ,
        bM ~ dnorm( 0 , 0.5 ) ,
        bA ~ dnorm( 0 , 0.5 ) ,
        sigma ~ dexp( 1 ),
      ## A -> M
        M ~ dnorm( mu_M , sigma_M ),
        mu_M <- aM + bAM*A,
        aM ~ dnorm( 0 , 0.2 ),
        bAM ~ dnorm( 0 , 0.5 ),
        sigma_M ~ dexp( 1 )
    ) , data = d )

test_that("R code 5.19",
    expect_equiv_eps( round(coef(m5.3_A),2) , 
        c( 0.00 ,  -0.07 ,  -0.61  ,  0.79 ,  0.00  , -0.69  ,  0.68 ) )
)

## R code 5.20
A_seq <- seq( from=-2 , to=2 , length.out=30 )

## R code 5.21
# prep data
sim_dat <- data.frame( A=A_seq )

# simulate M and then D, using A_seq
s <- sim( m5.3_A , data=sim_dat , vars=c("M","D") )

test_that("R code 5.21", {
    expect_named( s , c("M","D") )
    expect_equivalent( nrow(s$M) , 1000 )
    expect_equivalent( ncol(s$M) , 30 )
})

## R code 5.23
sim_dat <- data.frame( M=seq(from=-2,to=2,length.out=30) , A=0 )
s <- sim( m5.3_A , data=sim_dat , vars="D" )

test_that("R code 5.23", {
    expect_equivalent( dim(s) , c(1000,30) )
})

## R code 5.27
library(rethinking)
data(milk)
d <- milk

## R code 5.28
d$K <- scale( d$kcal.per.g )
d$N <- scale( d$neocortex.perc )
d$M <- scale( log(d$mass) )

## R code 5.29
test_that("R code 5.29 - NA test", {
expect_error( 
m5.5_draft <- quap(
    alist(
        K ~ dnorm( mu , sigma ) ,
        mu <- a + bN*N ,
        a ~ dnorm( 0 , 1 ) ,
        bN ~ dnorm( 0 , 1 ) ,
        sigma ~ dexp( 1 )
    ) , data=d )
) })

## R code 5.31
dcc <- d[ complete.cases(d$K,d$N,d$M) , ]

## R code 5.32

test_that("R code 5.32", {

    m5.5_draft <- quap(
    alist(
        K ~ dnorm( mu , sigma ) ,
        mu <- a + bN*N ,
        a ~ dnorm( 0 , 1 ) ,
        bN ~ dnorm( 0 , 1 ) ,
        sigma ~ dexp( 1 )
    ) , data=dcc )

    expect_equiv_eps( round(coef(m5.5_draft),2) , c( 0.09 , 0.16 , 1.00 ) )
})

## R code 5.34

test_that("R code 5.34", {

    m5.5 <- quap(
    alist(
        K ~ dnorm( mu , sigma ) ,
        mu <- a + bN*N ,
        a ~ dnorm( 0 , 0.2 ) ,
        bN ~ dnorm( 0 , 0.5 ) ,
        sigma ~ dexp( 1 )
    ) , data=dcc )

    expect_equiv_eps( round(coef(m5.5),2) , c(  0.04 , 0.13 , 1.00  ) )
})

## R code 5.37


test_that("R code 5.37", {

    m5.6 <- quap(
    alist(
        K ~ dnorm( mu , sigma ) ,
        mu <- a + bM*M ,
        a ~ dnorm( 0 , 0.2 ) ,
        bM ~ dnorm( 0 , 0.5 ) ,
        sigma ~ dexp( 1 )
    ) , data=dcc )

    expect_equiv_eps( round(coef(m5.6),2) , c( 0.05, -0.28,  0.95 ) )
})

## R code 5.38
test_that("R code 5.38", {

    m5.7 <- quap(
    alist(
        K ~ dnorm( mu , sigma ) ,
        mu <- a + bN*N + bM*M ,
        a ~ dnorm( 0 , 0.2 ) ,
        bN ~ dnorm( 0 , 0.5 ) ,
        bM ~ dnorm( 0 , 0.5 ) ,
        sigma ~ dexp( 1 )
    ) , data=dcc )

    expect_equiv_eps( round(coef(m5.7),2) , c( 0.07 , 0.68 ,-0.70 , 0.74 ) )
})

## R code 5.44
data(Howell1)
d <- Howell1

## R code 5.46
d$sex <- ifelse( d$male==1 , 2 , 1 )

## R code 5.47
test_that("R code 5.47", {

    m5.8 <- quap(
    alist(
        height ~ dnorm( mu , sigma ) ,
        mu <- a[sex] ,
        a[sex] ~ dnorm( 178 , 20 ) ,
        sigma ~ dunif( 0 , 50 )
    ) , data=d )

    expect_equiv_eps( round(coef(m5.8),2) , c(134.91, 142.58 , 27.31 ) )
})

## R code 5.49
data(milk)
d <- milk

## R code 5.50
d$clade_id <- as.integer( d$clade )

## R code 5.51
d$K <- scale( d$kcal.per.g )
test_that("R code 5.51", {

    m5.9 <- quap(
    alist(
        K ~ dnorm( mu , sigma ),
        mu <- a[clade_id],
        a[clade_id] ~ dnorm( 0 , 0.5 ),
        sigma ~ dexp( 1 )
    ) , data=d )

    expect_equiv_eps( round(coef(m5.9),2) , c(-0.48 , 0.37,  0.68, -0.59 , 0.72) )
})

## R code 5.52
set.seed(63)
d$house <- sample( rep(1:4,each=8) , size=nrow(d) )

## R code 5.53
test_that("R code 5.53", {

    set.seed(63)
    m5.10 <- quap(
    alist(
        K ~ dnorm( mu , sigma ),
        mu <- a[clade_id] + h[house],
        a[clade_id] ~ dnorm( 0 , 0.5 ),
        h[house] ~ dnorm( 0 , 0.5 ),
        sigma ~ dexp( 1 )
    ) , data=d )

    expect_equiv_eps( round(coef(m5.10),2) , 
        c(-0.42,  0.38 , 0.57, -0.51, -0.10 ,-0.20, -0.16 , 0.49 , 0.66) )
})

