# PSIS and WAIC tests
# test_dir('rethinking/tests/rethinking_tests',filter="PSIS_WAIC_compare")

context("PSIS, WAIC, and compare")

expect_equiv_eps <- function( x , y , eps=0.01 ) {
    expect_equivalent( x , y , tolerance=eps )
}

library(rethinking)
data(WaffleDivorce)
d <- WaffleDivorce

# standardize variables
d$A <- standardize( d$MedianAgeMarriage )
d$D <- standardize( d$Divorce )
d$M <- standardize( d$Marriage )

dat <- list( A=d$A , D=d$D , M=d$M )

set.seed(900)
m_quap <- quap(
    alist(
        D ~ dnorm( mu , sigma ) ,
        mu <- a + bM*M + bA*A ,
        a ~ dnorm( 0 , 0.2 ) ,
        bM ~ dnorm( 0 , 0.5 ) ,
        bA ~ dnorm( 0 , 0.5 ) ,
        sigma ~ dexp( 1 )
    ) , data = dat )

set.seed(900)
m_m2s <- map2stan(
    alist(
        D ~ dnorm( mu , sigma ) ,
        mu <- a + bM*M + bA*A ,
        a ~ dnorm( 0 , 0.2 ) ,
        bM ~ dnorm( 0 , 0.5 ) ,
        bA ~ dnorm( 0 , 0.5 ) ,
        sigma ~ dexp( 1 )
    ) , data = dat )

set.seed(900)
m_ulam <- ulam(
    alist(
        D ~ dnorm( mu , sigma ) ,
        mu <- a + bM*M + bA*A ,
        a ~ dnorm( 0 , 0.2 ) ,
        bM ~ dnorm( 0 , 0.5 ) ,
        bA ~ dnorm( 0 , 0.5 ) ,
        sigma ~ dexp( 1 )
    ) , data = dat , log_lik=TRUE )

post <- extract(m_ulam@stanfit)

test_that("PSIS quap",{
    set.seed(900)
    xxx <- PSIS(m_quap,pointwise=FALSE)
    expect_equiv_eps( round(xxx,2) , c(130.9 ,-65.45  ,  6.65  ,  15.8) )
    set.seed(900)
    xxx <- PSIS(m_quap,pointwise=TRUE)
    expect_equivalent( nrow(xxx) , 50 )
    expect_equivalent( ncol(xxx) , 5 )
    expect_equiv_eps( sum(xxx[[1]]) , 130.9 )
})

test_that("PSIS map2stan",{
    set.seed(900)
    xxx <- PSIS(m_m2s,pointwise=FALSE)
    expect_equiv_eps( round(xxx,2) , c(127.61 ,-63.81 ,   4.78  ,  12.9) )
    set.seed(900)
    xxx <- PSIS(m_m2s,pointwise=TRUE)
    expect_equivalent( nrow(xxx) , 50 )
    expect_equivalent( ncol(xxx) , 5 )
    expect_equiv_eps( sum(xxx[[1]]) , 127.6 )
})

test_that("PSIS ulam",{
    set.seed(900)
    xxx <- PSIS(m_ulam,pointwise=FALSE)
    expect_equiv_eps( round(xxx,2) , c(127.33 ,-63.66   , 4.55  , 12.58) )
    set.seed(900)
    xxx <- PSIS(m_ulam,pointwise=TRUE)
    expect_equivalent( nrow(xxx) , 50 )
    expect_equivalent( ncol(xxx) , 5 )
    expect_equiv_eps( sum(xxx[[1]]) , 128.33 )
})

test_that("PSIS list",{
    set.seed(900)
    xxx <- PSIS(post,pointwise=FALSE)
    expect_equiv_eps( round(xxx,2) , c(127.3 ,-63.65 ,   4.54 ,  12.57) )
    set.seed(900)
    xxx <- PSIS(post,pointwise=TRUE)
    expect_equivalent( nrow(xxx) , 50 )
    expect_equivalent( ncol(xxx) , 5 )
    expect_equiv_eps( sum(xxx[[1]]) , 127.3 )
})

test_that("WAIC quap",{
    set.seed(900)
    xxx <- WAIC(m_quap,pointwise=FALSE)
    expect_equiv_eps( round(xxx,2) , c(129.58 ,-58.8  ,  5.99  , 14.82) )
    set.seed(900)
    xxx <- WAIC(m_quap,pointwise=TRUE)
    expect_equivalent( nrow(xxx) , 50 )
    expect_equivalent( ncol(xxx) , 4 )
    expect_equiv_eps( sum(xxx[[1]]) , 129.58 )
})

test_that("WAIC map2stan",{
    set.seed(900)
    xxx <- WAIC(m_m2s,pointwise=FALSE)
    expect_equiv_eps( round(xxx,2) , c(127.43 ,-59.06  ,  4.66  , 12.65) )
    set.seed(900)
    xxx <- WAIC(m_m2s,pointwise=TRUE)
    expect_equivalent( nrow(xxx) , 50 )
    expect_equivalent( ncol(xxx) , 4 )
    expect_equiv_eps( sum(xxx[[1]]) , 127.4281 )
})

test_that("WAIC ulam",{
    set.seed(900)
    xxx <- WAIC(m_ulam,pointwise=FALSE)
    expect_equiv_eps( round(xxx,2) , c(127.21 ,-59.11  ,  4.49  , 12.41) )
    set.seed(900)
    xxx <- WAIC(m_ulam,pointwise=TRUE)
    expect_equivalent( nrow(xxx) , 50 )
    expect_equivalent( ncol(xxx) , 4 )
    expect_equiv_eps( sum(xxx[[1]]) , 127.21 )
})

test_that("WAIC list",{
    set.seed(900)
    xxx <- WAIC(post,pointwise=FALSE)
    expect_equiv_eps( round(xxx,2) , c(127.21, -59.11  ,  4.49  , 12.41) )
    set.seed(900)
    xxx <- WAIC(post,pointwise=TRUE)
    expect_equivalent( nrow(xxx) , 50 )
    expect_equivalent( ncol(xxx) , 4 )
    expect_equiv_eps( sum(xxx[[1]]) , 127.21 )
})

# compare()

set.seed(900)
m_quapM <- quap(
    alist(
        D ~ dnorm( mu , sigma ) ,
        mu <- a + bM*M ,
        a ~ dnorm( 0 , 0.2 ) ,
        bM ~ dnorm( 0 , 0.5 ) ,
        sigma ~ dexp( 1 )
    ) , data = dat )

set.seed(900)
m_quapA <- quap(
    alist(
        D ~ dnorm( mu , sigma ) ,
        mu <- a + bA*A ,
        a ~ dnorm( 0 , 0.2 ) ,
        bA ~ dnorm( 0 , 0.5 ) ,
        sigma ~ dexp( 1 )
    ) , data = dat )

test_that("compare",{
    set.seed(9000)
    xxx <- rethinking::compare( m_quap , m_quapM , m_quapA , func=WAIC )
    expect_equiv_eps( dim(xxx) , c(3,6) )
    xxx <- rethinking::compare( m_quap , m_quapM , m_quapA , func=PSIS )
    expect_equiv_eps( dim(xxx) , c(3,6) )
})


