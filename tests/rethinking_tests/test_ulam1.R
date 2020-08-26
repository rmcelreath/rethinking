# test_dir('rethinking/tests/rethinking_tests',filter="ulam1")

context('ulam')
library(rethinking)
# library(openssl) # for md5 hash

expect_equiv_eps <- function( x , y , eps=0.01 ) {
    expect_equivalent( x , y , tolerance=eps )
}

set_ulam_cmdstan(TRUE)

# student_t

data(Howell1)
d <- Howell1
d2 <- d[ d$age >= 18 , ]

set.seed(1)
mt <- ulam(
  alist(
    height ~ student_t( nu, mu , sigma ) ,
    mu <- a + b*weight ,
    a ~ normal( 178 , 20 ) ,
    b ~ normal( 0 , 1 ) ,
    sigma ~ exponential( 1 ),
    nu ~ gamma(2, 0.1)
  ) , data=d2 , sample=TRUE , seed=1 )

test_that("Student_t regression code match",
    expect_known_hash( mt@model , "98532fc851" )
)

coef_expect <- c(113.62  , 0.91  , 3.96  , 3.42  )
names(coef_expect) <- c("a","b","sigma","nu")
test_that("Student_t regression post means",
    expect_equiv_eps( round(coef(mt),2) , coef_expect )
)

rm(mt)
rm(d)
rm(d2)
gc()

# pareto

set.seed(49283)
y <- rpareto(1e3,1,1)
mp <- ulam(
    alist(
        y ~ pareto(1,alpha),
        alpha ~ exponential(1)
    ) , data=list(y=y) , chains=1 )

test_that("Pareto regression code match",
    expect_known_hash( mp@model , "7d9808c554" )
)

test_that("Pareto regression post match",
    expect_equiv_eps( round(coef(mp),2) , 1.01 )
)

# beta

test_that("ulam dbeta2 code",{
    y <- rbeta2( 100 , 0.3 , 3 )
    x <- rnorm( 100 )
    z <- ulam(
        alist(
            y ~ dbeta2( p , theta ),
            logit(p) <- a + b*x,
            a ~ normal(0,1.5),
            b ~ normal(0,0.5),
            theta ~ exponential(1)
        ), data=list(y=y,x=x) , sample=FALSE )
    expect_known_hash( z$model , "f2210bdc47" )
})

test_that("ulam dbeta2 no prior error",{
    y <- rbeta2( 100 , 0.3 , 3 )
    expect_error(
        z <- ulam(
        alist(
            y ~ dbeta2( p , theta ),
            logit(p) <- a,
            #a ~ normal(0,1.5),
            theta ~ exponential(1)
        ), data=list(y=y) , sample=TRUE )
    , "Do all parameters have priors?" )
})

# binomial tests

data( UCBadmit )
UCBadmit$male <- as.integer( ifelse( UCBadmit$applicant.gender=="male" , 1 , 0 ) )
UCBadmit$dept <- rep( 1:6 , each=2 )
UCBadmit$ag <- UCBadmit$applicant.gender
UCBadmit$applicant.gender <- NULL
UCBadmit$reject <- NULL

test_that("quap to ulam binomial",{
    z <- quap(
        alist(
            admit ~ dbinom(applications,p),
            logit(p) <- a,
            a ~ dnorm(0,4)
        ), data=UCBadmit )
    zz <- ulam( z , sample=FALSE )
    expect_known_hash( zz$model , "e1f6f2c9a3" )
})

test_that("ulam constraints",{
    z <- ulam(
        alist(
            admit ~ binomial(applications,p),
            logit(p) <- a,
            a ~ normal(0,4)
        ),
        constraints=list(a="lower=0"),
        data=UCBadmit , sample=FALSE )
    expect_known_hash( z$model , "608f2ad66b" )
})

test_that("ulam binomial log_lik",{
    z <- ulam(
        alist(
            admit ~ binomial(applications,p),
            logit(p) <- a + b*male,
            a ~ normal(0,1),
            b ~ normal(0,0.5)
        ), data=UCBadmit , sample=TRUE , log_lik=TRUE )
    expect_known_hash( z@model , "2add972dae" )
    expect_equivalent( length(WAIC(z)) , 4 )
})

# test to make sure constant '1' in binomial here isn't [i]'d in log_lik code
data( chimpanzees )
test_that("ulam binomial size 1 log_lik test",{
    z <- ulam(
        alist(
            y ~ binomial(1,p),
            logit(p) <- a,
            a ~ normal(0,4)
        ), data=list(y = chimpanzees$pulled_left) , sample=FALSE , log_lik=TRUE )
    expect_known_hash( z$model , "1eb0e4e07f" )
})

# continuous missing data

UCBadmit$x <- rnorm(12)
UCBadmit$x[1:2] <- NA
UCBadmit$x[2] <- NA

# explicit code with merge_missing
test_that("ulam continuous missing data 1",{
    z2b2 <- ulam(
        alist(
            admit ~ binomial(applications,p),
            logit(p) <- a + b*male + bx*x_merge,
            x_merge ~ normal( 0 , 1 ),
            x_merge <- merge_missing( x , x_impute ),
            a ~ normal(0,4),
            b ~ normal(0,1),
            bx ~ normal(0,1)
        ), data=UCBadmit , sample=TRUE , log_lik=FALSE )
    expect_known_hash( z2b2@model , "6300da1946" )
    expect_equivalent( dim(precis(z2b2,2)) , c(5,6) )
})

# automated version
# needs to build the merge code
test_that("ulam continuous missing data 2",{
    z2b_auto <- ulam(
        alist(
            admit ~ binomial(applications,p),
            logit(p) <- a + b*male + bx*x,
            x ~ normal( 0 , 1 ),
            a ~ normal(0,4),
            b ~ normal(0,1),
            bx ~ normal(0,1)
        ),
        data=UCBadmit , sample=TRUE , log_lik=TRUE )
    expect_known_hash( z2b_auto@model , "12312ebec6" )
    expect_equivalent( dim(precis(z2b_auto,2)) , c(5,6) )
})

UCBadmit$x2 <- rnorm(12)
UCBadmit$x2[5:7] <- NA

test_that("ulam continuous missing data 3",{
    z2c <- ulam(
        alist(
            admit ~ binomial(applications,p),
            logit(p) <- a + b*male + bx*x_merge + bx2*x2_merge,
            x_merge ~ normal( 0 , 1 ),
            x2_merge ~ normal( 0 , 1 ),
            x_merge <- merge_missing( x , x_impute ),
            x2_merge <- merge_missing( x2 , x2_impute ),
            a ~ normal(0,4),
            b ~ normal(0,1),
            c(bx,bx2) ~ normal(0,1)
        ), data=UCBadmit , sample=FALSE )
    expect_known_hash( z2c$model , "111badfdd3" )
})


# discrete missing data - need to mix over different terms
rm(UCBadmit)
data( UCBadmit )
UCBadmit$male <- ifelse( UCBadmit$applicant.gender=="male" , 1L , 0L )
UCBadmit$dept <- rep( 1:6 , each=2 )
UCBadmit$applicant.gender <- NULL
UCBadmit$male2 <- UCBadmit$male
UCBadmit$male2[1:2] <- (-9) # missingness code
UCBadmit$male2 <- as.integer(UCBadmit$male2)
UCBadmit$reject <- NULL

test_that("ulam discrete missing 1",{
    z <- ulam(
        alist(
            admit|male2==-9 ~ mixture(  
                phi_male , 
                binomial( applications , p_m1 ) , 
                binomial( applications , p_m0 ) ),
            admit|male2>-9 ~ binomial( applications , p ),
            logit(p) <- a[dept] + b*male2,
            logit(p_m1) <- a[dept] + b*1,
            logit(p_m0) <- a[dept] + b*0,
            male2|male2>-9 ~ bernoulli( phi_male ),
            phi_male ~ beta(2,2),
            a[dept] ~ normal(0,1),
            b ~ normal(0,1)
        ), data=UCBadmit , sample=TRUE , log_lik=FALSE )
    expect_known_hash( z@model , "aed7571a87" )
})

# same but with custom()
test_that("ulam discrete missing 2",{
    z2d <- ulam(
        alist(
            admit|male2==-1 ~ custom( log_mix( 
                phi_male , 
                binomial_lpmf(admit|applications,p_m1) , 
                binomial_lpmf(admit|applications,p_m0) ) ),
            admit|male2>-1 ~ binomial( applications , p ),
            logit(p) <- a[dept] + b*male2,
            logit(p_m1) <- a[dept] + b*1,
            logit(p_m0) <- a[dept] + b*0,
            male2|male2>-1 ~ bernoulli( phi_male ),
            phi_male ~ beta(2,2),
            a[dept] ~ normal(0,4),
            b ~ normal(0,1)
        ), data=UCBadmit , sample=FALSE , log_lik=FALSE )
    expect_known_hash( z2d$model , "73cd87c7f2" )
})

# log_sum_exp form
test_that("ulam discrete missing 3",{
    z2e <- ulam(
        alist(
            admit|male2==-1 ~ custom( 
                log_sum_exp( 
                    log(phi_male) + binomial_lpmf(admit|applications,p_m1) , 
                    log1m(phi_male) + binomial_lpmf(admit|applications,p_m0) ) ),
            admit|male2>-1 ~ binomial( applications , p ),
            logit(p) <- a + b*male2,
            logit(p_m1) <- a + b*1,
            logit(p_m0) <- a + b*0,
            male2|male2>-1 ~ bernoulli( phi_male ),
            phi_male ~ beta(2,2),
            a ~ normal(0,4),
            b ~ normal(0,1)
        ), data=UCBadmit , sample=FALSE )
    expect_known_hash( z2e$model , "dfca821d69" )
})

