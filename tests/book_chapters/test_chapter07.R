# test_dir('rethinking/tests/book_chapters',reporter="summary")
# test_dir('rethinking/tests/book_chapters',filter="chapter07")

context('chapter 7')
library(rethinking)

expect_equiv_eps <- function( x , y , eps=0.01 ) {
    expect_equivalent( x , y , tolerance=eps )
}


## R code 7.1
sppnames <- c( "afarensis","africanus","habilis","boisei",
    "rudolfensis","ergaster","sapiens")
brainvolcc <- c( 438 , 452 , 612, 521, 752, 871, 1350 )
masskg <- c( 37.0 , 35.5 , 34.5 , 41.5 , 55.5 , 61.0 , 53.5 )
d <- data.frame( species=sppnames , brain=brainvolcc , mass=masskg )

## R code 7.2
d$mass_std <- (d$mass - mean(d$mass))/sd(d$mass)
d$brain_std <- d$brain / max(d$brain)

## R code 7.3
set.seed(7.1)
m7.1 <- quap(
    alist(
        brain_std ~ dnorm( mu , exp(log_sigma) ),
        mu <- a + b*mass_std,
        a ~ dnorm( 0.5 , 1 ),
        b ~ dnorm( 0 , 10 ),
        log_sigma ~ dnorm( 0 , 1 )
    ), data=d )
test_that( "R code 7.3", {
    expect_equiv_eps( round(coef(m7.1),2) , c(0.53   ,   0.17  ,   -1.71) )
} )

## R code 7.6
R2_is_bad <- function( quap_fit ) {
    s <- sim( quap_fit , refresh=0 )
    r <- apply(s,2,mean) - d$brain_std
    1 - var2(r)/var2(d$brain_std)
}

## R code 7.7
m7.2 <- quap(
    alist(
        brain_std ~ dnorm( mu , exp(log_sigma) ),
        mu <- a + b[1]*mass_std + b[2]*mass_std^2,
        a ~ dnorm( 0.5 , 1 ),
        b ~ dnorm( 0 , 10 ),
        log_sigma ~ dnorm( 0 , 1 )
    ), data=d , start=list(b=rep(0,2)) )

## R code 7.8
m7.3 <- quap(
    alist(
        brain_std ~ dnorm( mu , exp(log_sigma) ),
        mu <- a + b[1]*mass_std + b[2]*mass_std^2 +
                  b[3]*mass_std^3,
        a ~ dnorm( 0.5 , 1 ),
        b ~ dnorm( 0 , 10 ),
        log_sigma ~ dnorm( 0 , 1 )
    ), data=d , start=list(b=rep(0,3)) )

m7.4 <- quap(
    alist(
        brain_std ~ dnorm( mu , exp(log_sigma) ),
        mu <- a + b[1]*mass_std + b[2]*mass_std^2 +
                  b[3]*mass_std^3 + b[4]*mass_std^4,
        a ~ dnorm( 0.5 , 1 ),
        b ~ dnorm( 0 , 10 ),
        log_sigma ~ dnorm( 0 , 1 )
    ), data=d , start=list(b=rep(0,4)) )

m7.5 <- quap(
    alist(
        brain_std ~ dnorm( mu , exp(log_sigma) ),
        mu <- a + b[1]*mass_std + b[2]*mass_std^2 +
                  b[3]*mass_std^3 + b[4]*mass_std^4 +
                  b[5]*mass_std^5,
        a ~ dnorm( 0.5 , 1 ),
        b ~ dnorm( 0 , 10 ),
        log_sigma ~ dnorm( 0 , 1 )
    ), data=d , start=list(b=rep(0,5)) )

## R code 7.9
m7.6 <- quap(
    alist(
        brain_std ~ dnorm( mu , 0.001 ),
        mu <- a + b[1]*mass_std + b[2]*mass_std^2 +
                  b[3]*mass_std^3 + b[4]*mass_std^4 +
                  b[5]*mass_std^5 + b[6]*mass_std^6,
        a ~ dnorm( 0.5 , 1 ),
        b ~ dnorm( 0 , 10 )
    ), data=d , start=list(b=rep(0,6)) )

## R code 7.10
post <- extract.samples(m7.1)
mass_seq <- seq( from=min(d$mass_std) , to=max(d$mass_std) , length.out=100 )
l <- link( m7.1 , data=list( mass_std=mass_seq ) )
test_that("R code 7.10", {
    expect_equivalent( dim(l) , c(1000,100) )
})

## R code 7.13

test_that("R code 7.13",{
    set.seed(1)
    xxx <- lppd( m7.1 , n=1e4 )
    expect_equivalent( length(xxx) , 7 )
    expect_equiv_eps( round(sum(xxx),1) , 2.5 )
})

## R code 7.14
test_that("R code 7.14",{
    set.seed(1)
    logprob <- sim( m7.1 , ll=TRUE , n=1e4 )
    expect_equiv_eps( round(sum(logprob),1) , 11672.9 )
})

## R code 7.15
test_that("R code 7.15",{
    set.seed(1)
    xxx <- sapply( list(m7.1,m7.2,m7.3,m7.4,m7.5,m7.6) , function(m) sum(lppd(m)) )
    expect_equiv_eps( sum(xxx) , 67.6713 )
})

## R code 6.13
set.seed(71)
# number of plants
N <- 100

# simulate initial heights
h0 <- rnorm(N,10,2)

# assign treatments and simulate fungus and growth
treatment <- rep( 0:1 , each=N/2 )
fungus <- rbinom( N , size=1 , prob=0.5 - treatment*0.4 )
h1 <- h0 + rnorm(N, 5 - 3*fungus)

# compose a clean data frame
d <- data.frame( h0=h0 , h1=h1 , treatment=treatment , fungus=fungus )

m6.6 <- quap(
    alist(
        h1 ~ dnorm( mu , sigma ),
        mu <- h0*p,
        p ~ dlnorm( 0 , 0.25 ),
        sigma ~ dexp( 1 )
    ), data=d )

m6.7 <- quap(
    alist(
        h1 ~ dnorm( mu , sigma ),
        mu <- h0 * p,
        p <- a + bt*treatment + bf*fungus,
        a ~ dlnorm( 0 , 0.2 ) ,
        bt ~ dnorm( 0 , 0.5 ),
        bf ~ dnorm( 0 , 0.5 ),
        sigma ~ dexp( 1 )
    ), data=d )

m6.8 <- quap(
    alist(
        h1 ~ dnorm( mu , sigma ),
        mu <- h0 * p,
        p <- a + bt*treatment,
        a ~ dlnorm( 0 , 0.2 ),
        bt ~ dnorm( 0 , 0.5 ),
        sigma ~ dexp( 1 )
    ), data=d )

## R code 7.25
test_that("R code 7.25",{
    set.seed(11)
    xxx <- WAIC( m6.7 )
    expect_equiv_eps( round(xxx,2) , c(361.45, -177.17,    3.55 ,  14.17) )
})

## R code 7.26
test_that("R code 7.26",{
    set.seed(77)
    xxx <- compare( m6.6 , m6.7 , m6.8 )
    expect_equiv_eps( dim(xxx) , c(3,6 ))
})

## R code 7.27
test_that("R code 7.27",{
    set.seed(91)
    waic_m6.7 <- WAIC( m6.7 , pointwise=TRUE )$WAIC
    waic_m6.8 <- WAIC( m6.8 , pointwise=TRUE )$WAIC
    n <- length(waic_m6.7)
    diff_m6.7_m6.8 <- waic_m6.7 - waic_m6.8
    xxx <- sqrt( n*var( diff_m6.7_m6.8 ) )
    expect_equiv_eps( xxx , 10.35785 )
})

## R code 7.30
test_that("R code 7.30",{
    set.seed(92)
    waic_m6.6 <- WAIC( m6.6 , pointwise=TRUE )$WAIC
    waic_m6.8 <- WAIC( m6.8 , pointwise=TRUE )$WAIC
    n <- length(waic_m6.6)
    diff_m6.6_m6.8 <- waic_m6.6 - waic_m6.8
    xxx <- sqrt( n*var( diff_m6.6_m6.8 ) )
    expect_equiv_eps( xxx , 4.858914 )
})

## R code 7.31
test_that("R code 7.31",{
    set.seed(93)
    xxx <- compare( m6.6 , m6.7 , m6.8 )@dSE
    expect_equiv_eps( dim(xxx) , c(3,3) )
})

## R code 7.32
library(rethinking)
data(WaffleDivorce)
d <- WaffleDivorce
d$A <- standardize( d$MedianAgeMarriage )
d$D <- standardize( d$Divorce )
d$M <- standardize( d$Marriage )

m5.1 <- quap(
    alist(
        D ~ dnorm( mu , sigma ) ,
        mu <- a + bA * A ,
        a ~ dnorm( 0 , 0.2 ) ,
        bA ~ dnorm( 0 , 0.5 ) ,
        sigma ~ dexp( 1 )
    ) , data = d )

m5.2 <- quap(
    alist(
        D ~ dnorm( mu , sigma ) ,
        mu <- a + bM * M ,
        a ~ dnorm( 0 , 0.2 ) ,
        bM ~ dnorm( 0 , 0.5 ) ,
        sigma ~ dexp( 1 )
    ) , data = d )

m5.3 <- quap(
    alist(
        D ~ dnorm( mu , sigma ) ,
        mu <- a + bM*M + bA*A ,
        a ~ dnorm( 0 , 0.2 ) ,
        bM ~ dnorm( 0 , 0.5 ) ,
        bA ~ dnorm( 0 , 0.5 ) ,
        sigma ~ dexp( 1 )
    ) , data = d )

## R code 7.33
test_that("R code 7.33",{
    set.seed(24071847)
    xxx <- compare( m5.1 , m5.2 , m5.3 , func=PSIS )
    expect_equiv_eps( dim(xxx) , c(3,6) )
})

## R code 7.34
test_that("R code 7.34",{
    set.seed(24071847)
    PSIS_m5.3 <- PSIS(m5.3,pointwise=TRUE)
    set.seed(24071847)
    WAIC_m5.3 <- WAIC(m5.3,pointwise=TRUE)
    expect_equiv_eps( dim(PSIS_m5.3) , c(50,5) )
    expect_equiv_eps( dim(WAIC_m5.3) , c(50,4) )
})

## R code 7.35
m5.3t <- quap(
    alist(
        D ~ dstudent( 2 , mu , sigma ) ,
        mu <- a + bM*M + bA*A ,
        a ~ dnorm( 0 , 0.2 ) ,
        bM ~ dnorm( 0 , 0.5 ) ,
        bA ~ dnorm( 0 , 0.5 ) ,
        sigma ~ dexp( 1 )
    ) , data = d )
