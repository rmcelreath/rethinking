# test_dir('rethinking/tests/book_chapters',reporter="summary")
# test_dir('rethinking/tests/book_chapters',filter="chapter15")

context('chapter 15')
library(rethinking)

expect_equiv_eps <- function( x , y , eps=0.01 ) {
    expect_equivalent( x , y , tolerance=eps )
}

## R code 15.2
library(rethinking)
data(WaffleDivorce)
d <- WaffleDivorce

## R code 15.3
dlist <- list(
    D_obs = standardize( d$Divorce ),
    D_sd = d$Divorce.SE / sd( d$Divorce ),
    M = standardize( d$Marriage ),
    A = standardize( d$MedianAgeMarriage ),
    N = nrow(d)
)

test_that("R code 15.3",{
    m15.1 <- ulam(
        alist(
            D_obs ~ dnorm( D_true , D_sd ),
            vector[N]:D_true ~ dnorm( mu , sigma ),
            mu <- a + bA*A + bM*M,
            a ~ dnorm(0,0.2),
            bA ~ dnorm(0,0.5),
            bM ~ dnorm(0,0.5),
            sigma ~ dexp(1)
        ) , data=dlist , chains=4 , cores=4 )
    expect_equivalent( dim(precis(m15.1,3)) , c(54,6) )
})

## R code 15.5
dlist <- list(
    D_obs = standardize( d$Divorce ),
    D_sd = d$Divorce.SE / sd( d$Divorce ),
    M_obs = standardize( d$Marriage ),
    M_sd = d$Marriage.SE / sd( d$Marriage ),
    A = standardize( d$MedianAgeMarriage ),
    N = nrow(d)
)

test_that("R code 15.5",{
    m15.2 <- ulam(
        alist(
            D_obs ~ dnorm( D_true , D_sd ),
            vector[N]:D_true ~ dnorm( mu , sigma ),
            mu <- a + bA*A + bM*M_true[i],
            M_obs ~ dnorm( M_true , M_sd ),
            vector[N]:M_true ~ dnorm( 0 , 1 ),
            a ~ dnorm(0,0.2),
            bA ~ dnorm(0,0.5),
            bM ~ dnorm(0,0.5),
            sigma ~ dexp( 1 )
        ) , data=dlist , chains=4 , cores=4 )
    expect_equivalent( dim(precis(m15.2,3)) , c(104,6) )
})

## R code 15.7
N <- 500
A <- rnorm(N)
M <- rnorm(N,-A)
D <- rnorm(N,A)
A_obs <- rnorm(N,A)

## R code 15.8
N <- 100
S <- rnorm( N )
H <- rbinom( N , size=10 , inv_logit(S) )

## R code 15.9
D <- rbern( N ) # dogs completely random
Hm <- H
Hm[D==1] <- NA

## R code 15.10
D <- ifelse( S > 0 , 1 , 0 )
Hm <- H
Hm[D==1] <- NA

## R code 15.11
set.seed(501)
N <- 1000
X <- rnorm(N)
S <- rnorm(N)
H <- rbinom( N , size=10 , inv_logit( 2 + S - 2*X ) )
D <- ifelse( X > 1 , 1 , 0 )
Hm <- H
Hm[D==1] <- NA

## R code 15.12
dat_list <- list(
    H = H,
    S = S )

test_that("R code 15.12",{
    m15.3 <- ulam(
        alist(
            H ~ binomial( 10 , p ),
            logit(p) <- a + bS*S,
            a ~ normal( 0 , 1 ),
            bS ~ normal( 0 , 0.5 )
        ), data=dat_list , chains=4 )
    expect_equivalent( dim(precis(m15.3,3)) , c(2,6) )
})

## R code 15.13
dat_list0 <- list(
    H = H[D==0],
    S = S[D==0] )

test_that("R code 15.13",{
    m15.4 <- ulam(
        alist(
            H ~ binomial( 10 , p ),
            logit(p) <- a + bS*S,
            a ~ normal( 0 , 1 ),
            bS ~ normal( 0 , 0.5 )
        ), data=dat_list0 , chains=4 )
    expect_equivalent( dim(precis(m15.4,3)) , c(2,6) )
})

## R code 15.16
library(rethinking)
data(milk)
d <- milk
d$neocortex.prop <- d$neocortex.perc / 100
d$logmass <- log(d$mass)
dat_list <- list(
    K = standardize( d$kcal.per.g ),
    B = standardize( d$neocortex.prop ),
    M = standardize( d$logmass ) )

## R code 15.17
test_that("R code 15.17",{
    m15.5 <- ulam(
        alist(
            K ~ dnorm( mu , sigma ),
            mu <- a + bB*B + bM*M,
            B ~ dnorm( nu , sigma_B ),
            c(a,nu) ~ dnorm( 0 , 0.5 ),
            c(bB,bM) ~ dnorm( 0, 0.5 ),
            sigma_B ~ dexp( 1 ),
            sigma ~ dexp( 1 )
        ) , data=dat_list , chains=4 , cores=4 )
    expect_equivalent( dim(precis(m15.5,3)) , c(18,6) )
})

## R code 15.19
obs_idx <- which( !is.na(d$neocortex.prop) )
dat_list_obs <- list(
    K = dat_list$K[obs_idx],
    B = dat_list$B[obs_idx],
    M = dat_list$M[obs_idx] )

test_that("R code 15.19",{
    m15.6 <- ulam(
        alist(
            K ~ dnorm( mu , sigma ),
            mu <- a + bB*B + bM*M,
            B ~ dnorm( nu , sigma_B ),
            c(a,nu) ~ dnorm( 0 , 0.5 ),
            c(bB,bM) ~ dnorm( 0, 0.5 ),
            sigma_B ~ dexp( 1 ),
            sigma ~ dexp( 1 )
        ) , data=dat_list_obs , chains=4 , cores=4 )
    expect_equivalent( dim(precis(m15.6,3)) , c(6,6) )
})

## R code 15.22
test_that("R code 15.22",{
    expect_warning(
    m15.7 <- ulam(
        alist(
           # K as function of B and M
            K ~ dnorm( mu , sigma ),
            mu <- a + bB*B_merge + bM*M,

           # M and B correlation
            MB ~ multi_normal( c(muM,muB) , Rho_BM , Sigma_BM ),
            matrix[29,2]:MB <<- append_col( M , B_merge ),

           # define B_merge as mix of observed and imputed values
            vector[29]:B_merge <- merge_missing( B , B_impute ),

           # priors
            c(a,muB,muM) ~ dnorm( 0 , 0.5 ),
            c(bB,bM) ~ dnorm( 0, 0.5 ),
            sigma ~ dexp( 1 ),
            Rho_BM ~ lkj_corr(2),
            Sigma_BM ~ dexp(1)
        ) , data=dat_list , chains=4 , cores=4 )
    )
    expect_equivalent( dim(precis(m15.7,3)) , c(24,6) )
})

## R code 15.29
set.seed(9)
N_houses <- 100L
alpha <- 5
beta <- (-3)
k <- 0.5
r <- 0.2
cat <- rbern( N_houses , k )
notes <- rpois( N_houses , alpha + beta*cat )
R_C <- rbern( N_houses , r )
cat_obs <- cat
cat_obs[R_C==1] <- (-9L)
dat <- list(
    notes = notes,
    cat = cat_obs,
    RC = R_C,
    N = as.integer(N_houses) )

## R code 15.30
test_that("R code 15.30",{
    m15.8 <- ulam(
        alist(
            # singing bird model
            ## cat known present/absent:
            notes|RC==0 ~ poisson( lambda ),
            log(lambda) <- a + b*cat,
            ## cat NA:
            notes|RC==1 ~ custom( log_sum_exp(
                    log(k) + poisson_lpmf( notes | exp(a + b) ),
                    log(1-k) + poisson_lpmf( notes | exp(a) )
                ) ),

            # priors
            a ~ normal(0,1),
            b ~ normal(0,0.5),

            # sneaking cat model
            cat|RC==0 ~ bernoulli(k),
            k ~ beta(2,2)
        ), data=dat , chains=4 , cores=4 )
    expect_equivalent( dim(precis(m15.8,3)) , c(3,6) )
})

## R code 15.31
test_that("R code 15.31",{
    m15.9 <- ulam(
        alist(
            # singing bird model
            notes|RC==0 ~ poisson( lambda ),
            notes|RC==1 ~ custom( log_sum_exp(
                    log(k) + poisson_lpmf( notes | exp(a + b) ),
                    log(1-k) + poisson_lpmf( notes | exp(a) )
                ) ),
            log(lambda) <- a + b*cat,
            a ~ normal(0,1),
            b ~ normal(0,0.5),

            # sneaking cat model
            cat|RC==0 ~ bernoulli(k),
            k ~ beta(2,2),

            # imputed values
            gq> vector[N]:PrC1 <- exp(lpC1)/(exp(lpC1)+exp(lpC0)),
            gq> vector[N]:lpC1 <- log(k) + poisson_lpmf( notes[i] | exp(a+b) ),
            gq> vector[N]:lpC0 <- log(1-k) + poisson_lpmf( notes[i] | exp(a) )
        ), data=dat , chains=4 , cores=4 )
    expect_equivalent( dim(precis(m15.9,3)) , c(303,6) )
})
