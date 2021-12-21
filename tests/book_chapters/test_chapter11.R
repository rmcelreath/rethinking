# test_dir('rethinking/tests/book_chapters',reporter="summary")
# test_dir('rethinking/tests/book_chapters',filter="chapter11")

context('chapter 11')
library(rethinking)

expect_equiv_eps <- function( x , y , eps=0.01 ) {
    expect_equivalent( x , y , tolerance=eps )
}


## R code 11.1
library(rethinking)
data(chimpanzees)
d <- chimpanzees

## R code 11.2
d$treatment <- 1 + d$prosoc_left + 2*d$condition

## R code 11.4
m11.1 <- quap(
    alist(
        pulled_left ~ dbinom( 1 , p ) ,
        logit(p) <- a ,
        a ~ dnorm( 0 , 10 )
    ) , data=d )

test_that("R code 11.4",{
    expect_equivalent( length(coef(m11.1)) , 1 )
})

## R code 11.5
set.seed(1999)
prior <- extract.prior( m11.1 , n=1e4 )

test_that("R code 11.5",{
    expect_equivalent( length(prior) , 1 )
    expect_equivalent( length(prior$a) , 10000 )
})

## R code 11.7
m11.2 <- quap(
    alist(
        pulled_left ~ dbinom( 1 , p ) ,
        logit(p) <- a + b[treatment] ,
        a ~ dnorm( 0 , 1.5 ),
        b[treatment] ~ dnorm( 0 , 10 )
    ) , data=d )
set.seed(1999)
prior <- extract.prior( m11.2 , n=1e4 )
p <- sapply( 1:4 , function(k) inv_logit( prior$a + prior$b[,k] ) )

test_that("R code 11.7",{
    expect_equivalent( length(prior) , 2 )
    expect_named( prior , c("a","b") )
    expect_equivalent( dim(prior$b) , c(10000,4) )
})

## R code 11.9
m11.3 <- quap(
    alist(
        pulled_left ~ dbinom( 1 , p ) ,
        logit(p) <- a + b[treatment] ,
        a ~ dnorm( 0 , 1.5 ),
        b[treatment] ~ dnorm( 0 , 0.5 )
    ) , data=d )
set.seed(1999)
prior <- extract.prior( m11.3 , n=1e4 )
p <- sapply( 1:4 , function(k) inv_logit( prior$a + prior$b[,k] ) )

test_that("R code 11.9",{
    expect_equivalent( length(prior) , 2 )
    expect_named( prior , c("a","b") )
    expect_equivalent( dim(prior$b) , c(10000,4) )
})

## R code 11.10
# prior trimmed data list
dat_list <- list(
    pulled_left = d$pulled_left,
    actor = d$actor,
    treatment = as.integer(d$treatment) )

## R code 11.11
m11.4 <- ulam(
    alist(
        pulled_left ~ dbinom( 1 , p ) ,
        logit(p) <- a[actor] + b[treatment] ,
        a[actor] ~ dnorm( 0 , 1.5 ),
        b[treatment] ~ dnorm( 0 , 0.5 )
    ) , data=dat_list , chains=4 , log_lik=TRUE )

test_that("R code 11.11",{
    xxx <- precis( m11.4 , depth=2 )
    expect_equivalent( dim(xxx) , c(11,6) )
})

## R code 11.12
post <- extract.samples(m11.4)

test_that("R code 11.12",{
    expect_equivalent( dim(post$a) , c(2000,7) )
    expect_equivalent( dim(post$b) , c(2000,4) )
})

## R code 11.18
d$side <- d$prosoc_left + 1 # right 1, left 2
d$cond <- d$condition + 1 # no partner 1, partner 2

## R code 11.19
dat_list2 <- list(
    pulled_left = d$pulled_left,
    actor = d$actor,
    side = d$side,
    cond = d$cond )
m11.5 <- ulam(
    alist(
        pulled_left ~ dbinom( 1 , p ) ,
        logit(p) <- a[actor] + bs[side] + bc[cond] ,
        a[actor] ~ dnorm( 0 , 1.5 ),
        bs[side] ~ dnorm( 0 , 0.5 ),
        bc[cond] ~ dnorm( 0 , 0.5 )
    ) , data=dat_list2 , chains=4 , log_lik=TRUE )

test_that("R code 11.19",{
    xxx <- precis( m11.5 , depth=2 )
    expect_equivalent( dim(xxx) , c(11,6) )
})

## R code 11.20
test_that("R code 11.20",{
    xxx <- compare( m11.5 , m11.4 , func=PSIS )
    expect_equivalent( dim(xxx) , c(2,6) )
})

## R code 11.21
post <- extract.samples( m11.4 , clean=FALSE )
test_that("R code 11.21",{
    expect_named( post , c('a', 'b', 'log_lik', 'p', 'lp__') )
    expect_equivalent( dim(post$log_lik) , c(2000,504) )
})

## R code 11.22
test_that("R code 11.22",{
    m11.4_stan_code <- stancode(m11.4)
    m11.4_stan <- cstan( model_code=m11.4_stan_code , data=dat_list , chains=4 )
    expect_warning( compare( m11.4_stan , m11.4 ) )
})

## R code 11.24
data(chimpanzees)
d <- chimpanzees
d$treatment <- 1 + d$prosoc_left + 2*d$condition
d$side <- d$prosoc_left + 1 # right 1, left 2
d$cond <- d$condition + 1 # no partner 1, partner 2
d_aggregated <- aggregate(
    d$pulled_left ,
    list( treatment=d$treatment , actor=d$actor ,
          side=d$side , cond=d$cond ) ,
    sum )
colnames(d_aggregated)[5] <- "left_pulls"

## R code 11.25
dat <- with( d_aggregated , list(
    left_pulls = left_pulls,
    treatment = treatment,
    actor = actor,
    side = side,
    cond = cond ) )

m11.6 <- ulam(
    alist(
        left_pulls ~ dbinom( 18 , p ) ,
        logit(p) <- a[actor] + b[treatment] ,
        a[actor] ~ dnorm( 0 , 1.5 ) ,
        b[treatment] ~ dnorm( 0 , 0.5 )
    ) , data=dat , chains=4 , log_lik=TRUE )

test_that("R code 11.25",{
    expect_equivalent( dim(precis(m11.6,2)) , c(11,6) )
})

## R code 11.26
test_that("R code 11.26",{
    expect_warning( rethinking::compare( m11.6 , m11.4 , func=PSIS ) )
})

## R code 11.28
library(rethinking)
data(UCBadmit)
d <- UCBadmit

## R code 11.29
dat_list <- list(
    admit = d$admit,
    applications = d$applications,
    gid = ifelse( d$applicant.gender=="male" , 1 , 2 )
)
m11.7 <- ulam(
    alist(
        admit ~ dbinom( applications , p ) ,
        logit(p) <- a[gid] ,
        a[gid] ~ dnorm( 0 , 1.5 )
    ) , data=dat_list , chains=4 )

test_that("R code 11.29",{
    expect_equivalent( dim(precis( m11.7 , depth=2 )) , c(2,6) )
})

## R code 11.31
test_that("R code 11.31",{
    xxx <- postcheck( m11.7 )
    expect_named( xxx , c("mean","PI","outPI") )
    expect_equivalent( dim(xxx$PI) , c(2,12) )
})

## R code 11.32
dat_list$dept_id <- rep(1:6,each=2)
m11.8 <- ulam(
    alist(
        admit ~ dbinom( applications , p ) ,
        logit(p) <- a[gid] + delta[dept_id] ,
        a[gid] ~ dnorm( 0 , 1.5 ) ,
        delta[dept_id] ~ dnorm( 0 , 1.5 )
    ) , data=dat_list , chains=4 , iter=4000 )

test_that("R code 11.32",{
    expect_equivalent( dim(precis( m11.8 , depth=2 )) , c(8,6) )
})

## R code 11.36
library(rethinking)
data(Kline)
d <- Kline

## R code 11.37
d$P <- scale( log(d$population) )
d$contact_id <- ifelse( d$contact=="high" , 2 , 1 )

## R code 11.45
dat <- list(
    T = d$total_tools ,
    P = d$P ,
    cid = d$contact_id )

# intercept only
m11.9 <- ulam(
    alist(
        T ~ dpois( lambda ),
        log(lambda) <- a,
        a ~ dnorm(3,0.5)
    ), data=dat , chains=4 , log_lik=TRUE )

test_that("R code 11.45",{
    expect_equivalent( dim(precis( m11.9 , depth=2 )) , c(1,6) )
})

# interaction model
m11.10 <- ulam(
    alist(
        T ~ dpois( lambda ),
        log(lambda) <- a[cid] + b[cid]*P,
        a[cid] ~ dnorm( 3 , 0.5 ),
        b[cid] ~ dnorm( 0 , 0.2 )
    ), data=dat , chains=4 , log_lik=TRUE )

test_that("R code 11.45",{
    expect_equivalent( dim(precis( m11.10 , depth=2 )) , c(4,6) )
})

## R code 11.47
test_that("R code 11.47",{
    expect_equivalent( length(PSIS( m11.10 , pointwise=TRUE )$k) , 10 )
})

## R code 11.49
dat2 <- list( T=d$total_tools, P=d$population, cid=d$contact_id )
m11.11 <- ulam(
    alist(
        T ~ dpois( lambda ),
        lambda <- exp(a[cid])*P^b[cid]/g,
        a[cid] ~ dnorm(1,1),
        b[cid] ~ dexp(1),
        g ~ dexp(1)
    ), data=dat2 , chains=4 , log_lik=TRUE )

test_that("R code 11.49",{
    expect_equivalent( dim(precis(m11.11,2)) , c(5,6) )
})

## R code 11.50
num_days <- 30
y <- rpois( num_days , 1.5 )

## R code 11.51
num_weeks <- 4
y_new <- rpois( num_weeks , 0.5*7 )

## R code 11.52
y_all <- c( y , y_new )
exposure <- c( rep(1,30) , rep(7,4) )
monastery <- c( rep(0,30) , rep(1,4) )
d <- data.frame( y=y_all , days=exposure , monastery=monastery )

## R code 11.53
# compute the offset
d$log_days <- log( d$days )

# fit the model
m11.12 <- quap(
    alist(
        y ~ dpois( lambda ),
        log(lambda) <- log_days + a + b*monastery,
        a ~ dnorm( 0 , 1 ),
        b ~ dnorm( 0 , 1 )
    ), data=d )

test_that("R code 11.53",{
    expect_equivalent( dim(precis(m11.12,2)) , c(2,4) )
})

## R code 11.54
post <- extract.samples( m11.12 )
lambda_old <- exp( post$a )
lambda_new <- exp( post$a + post$b )

test_that("R code 11.54",{
    xxx <- precis( data.frame( lambda_old , lambda_new ) )
    expect_equivalent( dim(xxx) , c(2,5) )
})

## R code 11.55
# simulate career choices among 500 individuals
N <- 500             # number of individuals
income <- c(1,2,5)   # expected income of each career
score <- 0.5*income  # scores for each career, based on income
# next line converts scores to probabilities
p <- softmax(score[1],score[2],score[3])

# now simulate choice
# outcome career holds event type values, not counts
career <- rep(NA,N)  # empty vector of choices for each individual
# sample chosen career for each individual
set.seed(34302)
for ( i in 1:N ) career[i] <- sample( 1:3 , size=1 , prob=p )

## R code 11.56
code_m11.13 <- "
data{
    int N; // number of individuals
    int K; // number of possible careers
    int career[N]; // outcome
    vector[K] career_income;
}
parameters{
    vector[K-1] a; // intercepts
    real<lower=0> b; // association of income with choice
}
model{
    vector[K] p;
    vector[K] s;
    a ~ normal( 0 , 1 );
    b ~ normal( 0 , 0.5 );
    s[1] = a[1] + b*career_income[1];
    s[2] = a[2] + b*career_income[2];
    s[3] = 0; // pivot
    p = softmax( s );
    career ~ categorical( p );
}
"

## R code 11.57
dat_list <- list( N=N , K=3 , career=career , career_income=income )

m11.13 <- cstan( model_code=code_m11.13 , data=dat_list , chains=4 )

test_that("R code 11.57",{
    
    expect_equivalent( dim( precis( m11.13 , 2 ) ) , c(3,6) )
})

## R code 11.58
post <- extract.samples( m11.13 )

# set up logit scores
s1 <- with( post , a[,1] + b*income[1] )
s2_orig <- with( post , a[,2] + b*income[2] )
s2_new <- with( post , a[,2] + b*income[2]*2 )

# compute probabilities for original and counterfactual
p_orig <- sapply( 1:length(post$b) , function(i)
    softmax( c(s1[i],s2_orig[i],0) ) )
p_new <- sapply( 1:length(post$b) , function(i)
    softmax( c(s1[i],s2_new[i],0) ) )

# summarize
p_diff <- p_new[2,] - p_orig[2,]

test_that("R code 11.58",{
    expect_equivalent( dim( precis( p_diff ) ) , c(1,5) )
})


## R code 11.59
N <- 500
# simulate family incomes for each individual
family_income <- runif(N)
# assign a unique coefficient for each type of event
b <- c(-2,0,2)
career <- rep(NA,N)  # empty vector of choices for each individual
for ( i in 1:N ) {
    score <- 0.5*(1:3) + b*family_income[i]
    p <- softmax(score[1],score[2],score[3])
    career[i] <- sample( 1:3 , size=1 , prob=p )
}

code_m11.14 <- "
data{
    int N; // number of observations
    int K; // number of outcome values
    int career[N]; // outcome
    real family_income[N];
}
parameters{
    vector[K-1] a; // intercepts
    vector[K-1] b; // coefficients on family income
}
model{
    vector[K] p;
    vector[K] s;
    a ~ normal(0,1.5);
    b ~ normal(0,1);
    for ( i in 1:N ) {
        for ( j in 1:(K-1) ) s[j] = a[j] + b[j]*family_income[i];
        s[K] = 0; // the pivot
        p = softmax( s );
        career[i] ~ categorical( p );
    }
}
"

dat_list <- list( N=N , K=3 , career=career , family_income=family_income )
m11.14 <- cstan( model_code=code_m11.14 , data=dat_list , chains=4 )

test_that("R code 11.59",{
    expect_equivalent( dim( precis( m11.14 , 2 ) ) , c(4,6) )
})


## R code 11.60
library(rethinking)
data(UCBadmit)
d <- UCBadmit

## R code 11.61
# binomial model of overall admission probability
m_binom <- quap(
    alist(
        admit ~ dbinom(applications,p),
        logit(p) <- a,
        a ~ dnorm( 0 , 1.5 )
    ), data=d )

# Poisson model of overall admission rate and rejection rate
# 'reject' is a reserved word in Stan, cannot use as variable name
dat <- list( admit=d$admit , rej=d$reject )
m_pois <- ulam(
    alist(
        admit ~ dpois(lambda1),
        rej ~ dpois(lambda2),
        log(lambda1) <- a1,
        log(lambda2) <- a2,
        c(a1,a2) ~ dnorm(0,1.5)
    ), data=dat , chains=3 , cores=3 )

## R code 11.66
library(rethinking)
data(AustinCats)
d <- AustinCats

d$adopt <- ifelse( d$out_event=="Adoption" , 1L , 0L )
dat <- list(
    days_to_event = as.numeric( d$days_to_event ),
    color_id = ifelse( d$color=="Black" , 1L , 2L ) ,
    adopted = d$adopt
)

m11.15 <- ulam(
    alist(
        days_to_event|adopted==1 ~ exponential( lambda ),
        days_to_event|adopted==0 ~ custom(exponential_lccdf( !Y | lambda )),
        lambda <- 1.0/mu,
        log(mu) <- a[color_id],
        a[color_id] ~ normal(0,1)
    ), data=dat , chains=4 , cores=4 )

test_that("R code 11.66",{
    expect_equivalent( dim( precis( m11.15 , 2 ) ) , c(2,6) )
})

## R code 11.67
post <- extract.samples( m11.15 )
post$D <- exp(post$a)
test_that("R code 11.67",{
    expect_equivalent( dim( precis( post , 2 ) ) , c(4,5) )
})


