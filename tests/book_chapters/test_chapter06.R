# test_dir('rethinking/tests/testthat',reporter="summary")
# test_dir('rethinking/tests/book_chapters',filter="chapter06")

context('chapter 6')
library(rethinking)
library(openssl) # for md5 hash

expect_equiv_eps <- function( x , y , eps=0.01 ) {
    expect_equivalent( x , y , tolerance=eps )
}

## R code 6.2
N <- 100                          # number of individuals
set.seed(909)
height <- rnorm(N,10,2)           # sim total height of each
leg_prop <- runif(N,0.4,0.5)      # leg as proportion of height
leg_left <- leg_prop*height +     # sim left leg as proportion + error
    rnorm( N , 0 , 0.02 )
leg_right <- leg_prop*height +    # sim right leg as proportion + error
    rnorm( N , 0 , 0.02 )
                                  # combine into data frame
d <- data.frame(height,leg_left,leg_right)

## R code 6.3
test_that("R code 6.3",{
    m6.1 <- quap(
    alist(
        height ~ dnorm( mu , sigma ) ,
        mu <- a + bl*leg_left + br*leg_right ,
        a ~ dnorm( 10 , 100 ) ,
        bl ~ dnorm( 2 , 10 ) ,
        br ~ dnorm( 2 , 10 ) ,
        sigma ~ dexp( 1 )
    ) ,
    data=d )
    expect_equiv_eps( round(coef(m6.1),2) , c( 0.98 , 0.21 , 1.78 , 0.62 ) )
})

## R code 6.7
test_that("R code 6.7",{
    set.seed(11001)
    m6.2 <- quap(
    alist(
        height ~ dnorm( mu , sigma ) ,
        mu <- a + bl*leg_left,
        a ~ dnorm( 10 , 100 ) ,
        bl ~ dnorm( 2 , 10 ) ,
        sigma ~ dexp( 1 )
    ) ,
    data=d )
    expect_equiv_eps( round(coef(m6.2),2) , c(  1.00 , 1.99 , 0.62 ) )
})

## R code 6.8
library(rethinking)
data(milk)
d <- milk
d$K <- scale( d$kcal.per.g )
d$F <- scale( d$perc.fat )
d$L <- scale( d$perc.lactose )

## R code 6.9
# kcal.per.g regressed on perc.fat
test_that("R code 6.9 m6.3",{
    m6.3 <- quap(
    alist(
        K ~ dnorm( mu , sigma ) ,
        mu <- a + bF*F ,
        a ~ dnorm( 0 , 0.2 ) ,
        bF ~ dnorm( 0 , 0.5 ) ,
        sigma ~ dexp( 1 )
    ) , data=d )
    expect_equiv_eps( round(coef(m6.3),2) , c( 0.00 , 0.86 , 0.45 ) )
})

# kcal.per.g regressed on perc.lactose
test_that("R code 6.9 m6.4",{
    m6.4 <- quap(
    alist(
        K ~ dnorm( mu , sigma ) ,
        mu <- a + bL*L ,
        a ~ dnorm( 0 , 0.2 ) ,
        bL ~ dnorm( 0 , 0.5 ) ,
        sigma ~ dexp( 1 )
    ) , data=d )
    expect_equiv_eps( round(coef(m6.4),2) , c( 0.00 ,-0.90 , 0.38  ) )
})

## R code 6.10
test_that("R code 6.10",{
    m6.5 <- quap(
    alist(
        K ~ dnorm( mu , sigma ) ,
        mu <- a + bF*F + bL*L ,
        a ~ dnorm( 0 , 0.2 ) ,
        bF ~ dnorm( 0 , 0.5 ) ,
        bL ~ dnorm( 0 , 0.5 ) ,
        sigma ~ dexp( 1 )
    ) ,
    data=d )
    expect_equiv_eps( round(coef(m6.5),2) , c( 0.00,  0.24 ,-0.68 , 0.38 ) )
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

## R code 6.15
test_that("R code 6.15",{
    m6.6 <- quap(
    alist(
        h1 ~ dnorm( mu , sigma ),
        mu <- h0*p,
        p ~ dlnorm( 0 , 0.25 ),
        sigma ~ dexp( 1 )
    ), data=d )
    expect_equiv_eps( round(coef(m6.6),2) , c( 1.43 , 1.79 ) )
})

## R code 6.16
test_that("R code 6.16",{
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
    expect_equiv_eps( round(coef(m6.7),2) , c( 1.48 , 0.00, -0.27 , 1.41 ) )
})

## R code 6.17
test_that("R code 6.17",{
    m6.8 <- quap(
    alist(
        h1 ~ dnorm( mu , sigma ),
        mu <- h0 * p,
        p <- a + bt*treatment,
        a ~ dlnorm( 0 , 0.2 ),
        bt ~ dnorm( 0 , 0.5 ),
        sigma ~ dexp( 1 )
    ), data=d )
    expect_equiv_eps( round(coef(m6.8),2) , c( 1.38 , 0.08 , 1.75 ) )
})

## R code 6.21
library(rethinking)
d <- sim_happiness( seed=1977 , N_years=1000 )

## R code 6.22
d2 <- d[ d$age>17 , ] # only adults
d2$A <- ( d2$age - 18 ) / ( 65 - 18 )

## R code 6.23
d2$mid <- d2$married + 1
test_that("R code 6.23",{
    m6.9 <- quap(
    alist(
        happiness ~ dnorm( mu , sigma ),
        mu <- a[mid] + bA*A,
        a[mid] ~ dnorm( 0 , 1 ),
        bA ~ dnorm( 0 , 2 ),
        sigma ~ dexp(1)
    ) , data=d2 )
    expect_equiv_eps( round(coef(m6.9),2) , c(-0.24 , 1.26 ,-0.75 , 0.99 ) )
})

## R code 6.24
test_that("R code 6.24",{
    m6.10 <- quap(
    alist(
        happiness ~ dnorm( mu , sigma ),
        mu <- a + bA*A,
        a ~ dnorm( 0 , 1 ),
        bA ~ dnorm( 0 , 2 ),
        sigma ~ dexp(1)
    ) , data=d2 )
    expect_equiv_eps( round(coef(m6.10),2) , c( 0.00 , 0.00 , 1.21 ) )
})

## R code 6.25
N <- 200  # number of grandparent-parent-child triads
b_GP <- 1 # direct effect of G on P
b_GC <- 0 # direct effect of G on C
b_PC <- 1 # direct effect of P on C
b_U <- 2  # direct effect of U on P and C

## R code 6.26
set.seed(1)
U <- 2*rbern( N , 0.5 ) - 1
G <- rnorm( N )
P <- rnorm( N , b_GP*G + b_U*U )
C <- rnorm( N , b_PC*P + b_GC*G + b_U*U )
d <- data.frame( C=C , P=P , G=G , U=U )

## R code 6.27
test_that("R code 6.27",{
    m6.11 <- quap(
    alist(
        C ~ dnorm( mu , sigma ),
        mu <- a + b_PC*P + b_GC*G,
        a ~ dnorm( 0 , 1 ),
        c(b_PC,b_GC) ~ dnorm( 0 , 1 ),
        sigma ~ dexp( 1 )
    ), data=d )
    expect_equiv_eps( round(coef(m6.11),2) , c(-0.12 , 1.79 ,-0.84 , 1.41 ) )
})

## R code 6.28
test_that("R code 6.28",{
    m6.12 <- quap(
    alist(
        C ~ dnorm( mu , sigma ),
        mu <- a + b_PC*P + b_GC*G + b_U*U,
        a ~ dnorm( 0 , 1 ),
        c(b_PC,b_GC,b_U) ~ dnorm( 0 , 1 ),
        sigma ~ dexp( 1 )
    ), data=d )
    expect_equiv_eps( round(coef(m6.12),2) , 
        c(-0.12 , 1.01 ,-0.04 , 2.00 , 1.02 ) )
})
