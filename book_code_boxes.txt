## R code 0.1
print( "All models are wrong, but some are useful." )

## R code 0.2
x <- 1:2
x <- x*10
x <- log(x)
x <- sum(x)
x <- exp(x)
x

## R code 0.3
( log( 0.01^200 ) )
( 200 * log(0.01) )

## R code 0.4
# Load the data:
# car braking distances in feet paired with speeds in km/h
# see ?cars for details
data(cars)

# fit a linear regression of distance on speed
m <- lm( dist ~ speed , data=cars )

# estimated coefficients from the model
coef(m)

# plot residuals against speed
plot( resid(m) ~ speed , data=cars )

## R code 0.5
install.packages(c("coda","mvtnorm","devtools"))
library(devtools)
devtools::install_github("rmcelreath/rethinking")

## R code 2.1
ways <- c( 0 , 3 , 8 , 9 , 0 )
ways/sum(ways)

## R code 2.2
dbinom( 6 , size=9 , prob=0.5 )

## R code 2.3
# define grid
p_grid <- seq( from=0 , to=1 , length.out=20 )

# define prior
prior <- rep( 1 , 20 )

# compute likelihood at each value in grid
likelihood <- dbinom( 6 , size=9 , prob=p_grid )

# compute product of likelihood and prior
unstd.posterior <- likelihood * prior

# standardize the posterior, so it sums to 1
posterior <- unstd.posterior / sum(unstd.posterior)

## R code 2.4
plot( p_grid , posterior , type="b" ,
    xlab="probability of water" , ylab="posterior probability" )
mtext( "20 points" )

## R code 2.5
prior <- ifelse( p_grid < 0.5 , 0 , 1 )
prior <- exp( -5*abs( p_grid - 0.5 ) )

## R code 2.6
library(rethinking)
globe.qa <- map(
    alist(
        w ~ dbinom(9,p) ,  # binomial likelihood
        p ~ dunif(0,1)     # uniform prior
    ) ,
    data=list(w=6) )

# display summary of quadratic approximation
precis( globe.qa )

## R code 2.7
# analytical calculation
w <- 6
n <- 9
curve( dbeta( x , w+1 , n-w+1 ) , from=0 , to=1 )
# quadratic approximation
curve( dnorm( x , 0.67 , 0.16 ) , lty=2 , add=TRUE )

## R code 3.1
PrPV <- 0.95
PrPM <- 0.01
PrV <- 0.001
PrP <- PrPV*PrV + PrPM*(1-PrV)
( PrVP <- PrPV*PrV / PrP )

## R code 3.2
p_grid <- seq( from=0 , to=1 , length.out=1000 )
prior <- rep( 1 , 1000 )
likelihood <- dbinom( 6 , size=9 , prob=p_grid )
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)

## R code 3.3
samples <- sample( p_grid , prob=posterior , size=1e4 , replace=TRUE )

## R code 3.4
plot( samples )

## R code 3.5
library(rethinking)
dens( samples )

## R code 3.6
# add up posterior probability where p < 0.5
sum( posterior[ p_grid < 0.5 ] )

## R code 3.7
sum( samples < 0.5 ) / 1e4

## R code 3.8
sum( samples > 0.5 & samples < 0.75 ) / 1e4

## R code 3.9
quantile( samples , 0.8 )

## R code 3.10
quantile( samples , c( 0.1 , 0.9 ) )

## R code 3.11
p_grid <- seq( from=0 , to=1 , length.out=1000 )
prior <- rep(1,1000)
likelihood <- dbinom( 3 , size=3 , prob=p_grid )
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)
samples <- sample( p_grid , size=1e4 , replace=TRUE , prob=posterior )

## R code 3.12
PI( samples , prob=0.5 )

## R code 3.13
HPDI( samples , prob=0.5 )

## R code 3.14
p_grid[ which.max(posterior) ]

## R code 3.15
chainmode( samples , adj=0.01 )

## R code 3.16
mean( samples )
median( samples )

## R code 3.17
sum( posterior*abs( 0.5 - p_grid ) )

## R code 3.18
loss <- sapply( p_grid , function(d) sum( posterior*abs( d - p_grid ) ) )

## R code 3.19
p_grid[ which.min(loss) ]

## R code 3.20
dbinom( 0:2 , size=2 , prob=0.7 )

## R code 3.21
rbinom( 1 , size=2 , prob=0.7 )

## R code 3.22
rbinom( 10 , size=2 , prob=0.7 )

## R code 3.23
dummy_w <- rbinom( 1e5 , size=2 , prob=0.7 )
table(dummy_w)/1e5

## R code 3.24
dummy_w <- rbinom( 1e5 , size=9 , prob=0.7 )
simplehist( dummy_w , xlab="dummy water count" )

## R code 3.25
w <- rbinom( 1e4 , size=9 , prob=0.6 )

## R code 3.26
w <- rbinom( 1e4 , size=9 , prob=samples )

## R code 3.27
p_grid <- seq( from=0 , to=1 , length.out=1000 )
prior <- rep( 1 , 1000 )
likelihood <- dbinom( 6 , size=9 , prob=p_grid )
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)
set.seed(100)
samples <- sample( p_grid , prob=posterior , size=1e4 , replace=TRUE )

## R code 3.28
birth1 <- c(1,0,0,0,1,1,0,1,0,1,0,0,1,1,0,1,1,0,0,0,1,0,0,0,1,0,
0,0,0,1,1,1,0,1,0,1,1,1,0,1,0,1,1,0,1,0,0,1,1,0,1,0,0,0,0,0,0,0,
1,1,0,1,0,0,1,0,0,0,1,0,0,1,1,1,1,0,1,0,1,1,1,1,1,0,0,1,0,1,1,0,
1,0,1,1,1,0,1,1,1,1)
birth2 <- c(0,1,0,1,0,1,1,1,0,0,1,1,1,1,1,0,0,1,1,1,0,0,1,1,1,0,
1,1,1,0,1,1,1,0,1,0,0,1,1,1,1,0,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,
1,1,1,0,1,1,0,1,1,0,1,1,1,0,0,0,0,0,0,1,0,0,0,1,1,0,0,1,0,0,1,1,
0,0,0,1,1,1,0,0,0,0)

## R code 3.29
library(rethinking)
data(homeworkch3)

## R code 3.30
sum(birth1) + sum(birth2)

## R code 4.1
pos <- replicate( 1000 , sum( runif(16,-1,1) ) )

## R code 4.2
prod( 1 + runif(12,0,0.1) )

## R code 4.3
growth <- replicate( 10000 , prod( 1 + runif(12,0,0.1) ) )
dens( growth , norm.comp=TRUE )

## R code 4.4
big <- replicate( 10000 , prod( 1 + runif(12,0,0.5) ) )
small <- replicate( 10000 , prod( 1 + runif(12,0,0.01) ) )

## R code 4.5
log.big <- replicate( 10000 , log(prod(1 + runif(12,0,0.5))) )

## R code 4.6
w <- 6; n <- 9;
p_grid <- seq(from=0,to=1,length.out=100)
posterior <- dbinom(w,n,p_grid)*dunif(p_grid,0,1)
posterior <- posterior/sum(posterior)

## R code 4.7
library(rethinking)
data(Howell1)
d <- Howell1

## R code 4.8
str( d )

## R code 4.9
d$height

## R code 4.10
d2 <- d[ d$age >= 18 , ]

## R code 4.11
curve( dnorm( x , 178 , 20 ) , from=100 , to=250 )

## R code 4.12
curve( dunif( x , 0 , 50 ) , from=-10 , to=60 )

## R code 4.13
sample_mu <- rnorm( 1e4 , 178 , 20 )
sample_sigma <- runif( 1e4 , 0 , 50 )
prior_h <- rnorm( 1e4 , sample_mu , sample_sigma )
dens( prior_h )

## R code 4.14
mu.list <- seq( from=140, to=160 , length.out=200 )
sigma.list <- seq( from=4 , to=9 , length.out=200 )
post <- expand.grid( mu=mu.list , sigma=sigma.list )
post$LL <- sapply( 1:nrow(post) , function(i) sum( dnorm(
                d2$height ,
                mean=post$mu[i] ,
                sd=post$sigma[i] ,
                log=TRUE ) ) )
post$prod <- post$LL + dnorm( post$mu , 178 , 20 , TRUE ) +
    dunif( post$sigma , 0 , 50 , TRUE )
post$prob <- exp( post$prod - max(post$prod) )

## R code 4.15
contour_xyz( post$mu , post$sigma , post$prob )

## R code 4.16
image_xyz( post$mu , post$sigma , post$prob )

## R code 4.17
sample.rows <- sample( 1:nrow(post) , size=1e4 , replace=TRUE ,
    prob=post$prob )
sample.mu <- post$mu[ sample.rows ]
sample.sigma <- post$sigma[ sample.rows ]

## R code 4.18
plot( sample.mu , sample.sigma , cex=0.5 , pch=16 , col=col.alpha(rangi2,0.1) )

## R code 4.19
dens( sample.mu )
dens( sample.sigma )

## R code 4.20
HPDI( sample.mu )
HPDI( sample.sigma )

## R code 4.21
d3 <- sample( d2$height , size=20 )

## R code 4.22
mu.list <- seq( from=150, to=170 , length.out=200 )
sigma.list <- seq( from=4 , to=20 , length.out=200 )
post2 <- expand.grid( mu=mu.list , sigma=sigma.list )
post2$LL <- sapply( 1:nrow(post2) , function(i)
    sum( dnorm( d3 , mean=post2$mu[i] , sd=post2$sigma[i] ,
    log=TRUE ) ) )
post2$prod <- post2$LL + dnorm( post2$mu , 178 , 20 , TRUE ) +
    dunif( post2$sigma , 0 , 50 , TRUE )
post2$prob <- exp( post2$prod - max(post2$prod) )
sample2.rows <- sample( 1:nrow(post2) , size=1e4 , replace=TRUE ,
    prob=post2$prob )
sample2.mu <- post2$mu[ sample2.rows ]
sample2.sigma <- post2$sigma[ sample2.rows ]
plot( sample2.mu , sample2.sigma , cex=0.5 ,
    col=col.alpha(rangi2,0.1) ,
    xlab="mu" , ylab="sigma" , pch=16 )

## R code 4.23
dens( sample2.sigma , norm.comp=TRUE )

## R code 4.24
library(rethinking)
data(Howell1)
d <- Howell1
d2 <- d[ d$age >= 18 , ]

## R code 4.25
flist <- alist(
    height ~ dnorm( mu , sigma ) ,
    mu ~ dnorm( 178 , 20 ) ,
    sigma ~ dunif( 0 , 50 )
)

## R code 4.26
m4.1 <- map( flist , data=d2 )

## R code 4.27
precis( m4.1 )

## R code 4.28
start <- list(
    mu=mean(d2$height),
    sigma=sd(d2$height)
)

## R code 4.29
m4.2 <- map(
        alist(
            height ~ dnorm( mu , sigma ) ,
            mu ~ dnorm( 178 , 0.1 ) ,
            sigma ~ dunif( 0 , 50 )
        ) ,
        data=d2 )
precis( m4.2 )

## R code 4.30
vcov( m4.1 )

## R code 4.31
diag( vcov( m4.1 ) )
cov2cor( vcov( m4.1 ) )

## R code 4.32
library(rethinking)
post <- extract.samples( m4.1 , n=1e4 )
head(post)

## R code 4.33
precis(post)

## R code 4.34
library(MASS)
post <- mvrnorm( n=1e4 , mu=coef(m4.1) , Sigma=vcov(m4.1) )

## R code 4.35
m4.1_logsigma <- map(
        alist(
            height ~ dnorm( mu , exp(log_sigma) ) ,
            mu ~ dnorm( 178 , 20 ) ,
            log_sigma ~ dnorm( 2 , 10 )
        ) , data=d2 )

## R code 4.36
post <- extract.samples( m4.1_logsigma )
sigma <- exp( post$log_sigma )

## R code 4.37
plot( d2$height ~ d2$weight )

## R code 4.38
# load data again, since it's a long way back
library(rethinking)
data(Howell1)
d <- Howell1
d2 <- d[ d$age >= 18 , ]

# fit model
m4.3 <- map(
    alist(
        height ~ dnorm( mu , sigma ) ,
        mu <- a + b*weight ,
        a ~ dnorm( 156 , 100 ) ,
        b ~ dnorm( 0 , 10 ) ,
        sigma ~ dunif( 0 , 50 )
    ) ,
    data=d2 )

## R code 4.39
m4.3 <- map(
    alist(
        height ~ dnorm( a + b*weight , sigma ) ,
        a ~ dnorm( 178 , 100 ) ,
        b ~ dnorm( 0 , 10 ) ,
        sigma ~ dunif( 0 , 50 )
    ) ,
    data=d2 )

## R code 4.40
precis( m4.3 )

## R code 4.41
precis( m4.3 , corr=TRUE )

## R code 4.42
d2$weight.c <- d2$weight - mean(d2$weight)

## R code 4.43
m4.4 <- map(
    alist(
        height ~ dnorm( mu , sigma ) ,
        mu <- a + b*weight.c ,
        a ~ dnorm( 178 , 100 ) ,
        b ~ dnorm( 0 , 10 ) ,
        sigma ~ dunif( 0 , 50 )
    ) ,
    data=d2 )

## R code 4.44
precis( m4.4 , corr=TRUE )

## R code 4.45
plot( height ~ weight , data=d2 )
abline( a=coef(m4.3)["a"] , b=coef(m4.3)["b"] )

## R code 4.46
post <- extract.samples( m4.3 )

## R code 4.47
post[1:5,]

## R code 4.48
N <- 10
dN <- d2[ 1:N , ]
mN <- map(
    alist(
        height ~ dnorm( mu , sigma ) ,
        mu <- a + b*weight ,
        a ~ dnorm( 178 , 100 ) ,
        b ~ dnorm( 0 , 10 ) ,
        sigma ~ dunif( 0 , 50 )
    ) , data=dN )

## R code 4.49
# extract 20 samples from the posterior
post <- extract.samples( mN , n=20 )

# display raw data and sample size
plot( dN$weight , dN$height ,
    xlim=range(d2$weight) , ylim=range(d2$height) ,
    col=rangi2 , xlab="weight" , ylab="height" )
mtext(concat("N = ",N))

# plot the lines, with transparency
for ( i in 1:20 )
    abline( a=post$a[i] , b=post$b[i] , col=col.alpha("black",0.3) )

## R code 4.50
mu_at_50 <- post$a + post$b * 50

## R code 4.51
dens( mu_at_50 , col=rangi2 , lwd=2 , xlab="mu|weight=50" )

## R code 4.52
HPDI( mu_at_50 , prob=0.89 )

## R code 4.53
mu <- link( m4.3 )
str(mu)

## R code 4.54
# define sequence of weights to compute predictions for
# these values will be on the horizontal axis
weight.seq <- seq( from=25 , to=70 , by=1 )

# use link to compute mu
# for each sample from posterior
# and for each weight in weight.seq
mu <- link( m4.3 , data=data.frame(weight=weight.seq) )
str(mu)

## R code 4.55
# use type="n" to hide raw data
plot( height ~ weight , d2 , type="n" )

# loop over samples and plot each mu value
for ( i in 1:100 )
    points( weight.seq , mu[i,] , pch=16 , col=col.alpha(rangi2,0.1) )

## R code 4.56
# summarize the distribution of mu
mu.mean <- apply( mu , 2 , mean )
mu.HPDI <- apply( mu , 2 , HPDI , prob=0.89 )

## R code 4.57
# plot raw data
# fading out points to make line and interval more visible
plot( height ~ weight , data=d2 , col=col.alpha(rangi2,0.5) )

# plot the MAP line, aka the mean mu for each weight
lines( weight.seq , mu.mean )

# plot a shaded region for 89% HPDI
shade( mu.HPDI , weight.seq )

## R code 4.58
post <- extract.samples(m4.3)
mu.link <- function(weight) post$a + post$b*weight
weight.seq <- seq( from=25 , to=70 , by=1 )
mu <- sapply( weight.seq , mu.link )
mu.mean <- apply( mu , 2 , mean )
mu.HPDI <- apply( mu , 2 , HPDI , prob=0.89 )

## R code 4.59
sim.height <- sim( m4.3 , data=list(weight=weight.seq) )
str(sim.height)

## R code 4.60
height.PI <- apply( sim.height , 2 , PI , prob=0.89 )

## R code 4.61
# plot raw data
plot( height ~ weight , d2 , col=col.alpha(rangi2,0.5) )

# draw MAP line
lines( weight.seq , mu.mean )

# draw HPDI region for line
shade( mu.HPDI , weight.seq )

# draw PI region for simulated heights
shade( height.PI , weight.seq )

## R code 4.62
sim.height <- sim( m4.3 , data=list(weight=weight.seq) , n=1e4 )
height.PI <- apply( sim.height , 2 , PI , prob=0.89 )

## R code 4.63
post <- extract.samples(m4.3)
weight.seq <- 25:70
sim.height <- sapply( weight.seq , function(weight)
    rnorm(
        n=nrow(post) ,
        mean=post$a + post$b*weight ,
        sd=post$sigma ) )
height.PI <- apply( sim.height , 2 , PI , prob=0.89 )

## R code 4.64
library(rethinking)
data(Howell1)
d <- Howell1
str(d)

## R code 4.65
d$weight.s <- ( d$weight - mean(d$weight) )/sd(d$weight)

## R code 4.66
d$weight.s2 <- d$weight.s^2
m4.5 <- map(
    alist(
        height ~ dnorm( mu , sigma ) ,
        mu <- a + b1*weight.s + b2*weight.s2 ,
        a ~ dnorm( 178 , 100 ) ,
        b1 ~ dnorm( 0 , 10 ) ,
        b2 ~ dnorm( 0 , 10 ) ,
        sigma ~ dunif( 0 , 50 )
    ) ,
    data=d )

## R code 4.67
precis( m4.5 )

## R code 4.68
weight.seq <- seq( from=-2.2 , to=2 , length.out=30 )
pred_dat <- list( weight.s=weight.seq , weight.s2=weight.seq^2 )
mu <- link( m4.5 , data=pred_dat )
mu.mean <- apply( mu , 2 , mean )
mu.PI <- apply( mu , 2 , PI , prob=0.89 )
sim.height <- sim( m4.5 , data=pred_dat )
height.PI <- apply( sim.height , 2 , PI , prob=0.89 )

## R code 4.69
plot( height ~ weight.s , d , col=col.alpha(rangi2,0.5) )
lines( weight.seq , mu.mean )
shade( mu.PI , weight.seq )
shade( height.PI , weight.seq )

## R code 4.70
d$weight.s3 <- d$weight.s^3
m4.6 <- map(
    alist(
        height ~ dnorm( mu , sigma ) ,
        mu <- a + b1*weight.s + b2*weight.s2 + b3*weight.s3 ,
        a ~ dnorm( 178 , 100 ) ,
        b1 ~ dnorm( 0 , 10 ) ,
        b2 ~ dnorm( 0 , 10 ) ,
        b3 ~ dnorm( 0 , 10 ) ,
        sigma ~ dunif( 0 , 50 )
    ) ,
    data=d )

## R code 4.71
plot( height ~ weight.s , d , col=col.alpha(rangi2,0.5) , xaxt="n" )

## R code 4.72
at <- c(-2,-1,0,1,2)
labels <- at*sd(d$weight) + mean(d$weight)
axis( side=1 , at=at , labels=round(labels,1) )

## R code 4.73
plot( height ~ weight , data=Howell1 ,
    col=col.alpha(rangi2,0.4) )

## R code 5.1
# load data
library(rethinking)
data(WaffleDivorce)
d <- WaffleDivorce

# standardize predictor
d$MedianAgeMarriage.s <- (d$MedianAgeMarriage-mean(d$MedianAgeMarriage))/
    sd(d$MedianAgeMarriage)

# fit model
m5.1 <- map(
    alist(
        Divorce ~ dnorm( mu , sigma ) ,
        mu <- a + bA * MedianAgeMarriage.s ,
        a ~ dnorm( 10 , 10 ) ,
        bA ~ dnorm( 0 , 1 ) ,
        sigma ~ dunif( 0 , 10 )
    ) , data = d )

## R code 5.2
# compute percentile interval of mean
MAM.seq <- seq( from=-3 , to=3.5 , length.out=30 )
mu <- link( m5.1 , data=data.frame(MedianAgeMarriage.s=MAM.seq) )
mu.PI <- apply( mu , 2 , PI )

# plot it all
plot( Divorce ~ MedianAgeMarriage.s , data=d , col=rangi2 )
abline( m5.1 )
shade( mu.PI , MAM.seq )

## R code 5.3
d$Marriage.s <- (d$Marriage - mean(d$Marriage))/sd(d$Marriage)
m5.2 <- map(
    alist(
        Divorce ~ dnorm( mu , sigma ) ,
        mu <- a + bR * Marriage.s ,
        a ~ dnorm( 10 , 10 ) ,
        bR ~ dnorm( 0 , 1 ) ,
        sigma ~ dunif( 0 , 10 )
    ) , data = d )

## R code 5.4
m5.3 <- map(
    alist(
        Divorce ~ dnorm( mu , sigma ) ,
        mu <- a + bR*Marriage.s + bA*MedianAgeMarriage.s ,
        a ~ dnorm( 10 , 10 ) ,
        bR ~ dnorm( 0 , 1 ) ,
        bA ~ dnorm( 0 , 1 ) ,
        sigma ~ dunif( 0 , 10 )
    ) ,
    data = d )
precis( m5.3 )

## R code 5.5
plot( precis(m5.3) )

## R code 5.6
m5.4 <- map(
    alist(
        Marriage.s ~ dnorm( mu , sigma ) ,
        mu <- a + b*MedianAgeMarriage.s ,
        a ~ dnorm( 0 , 10 ) ,
        b ~ dnorm( 0 , 1 ) ,
        sigma ~ dunif( 0 , 10 )
    ) ,
    data = d )

## R code 5.7
# compute expected value at MAP, for each State
mu <- coef(m5.4)['a'] + coef(m5.4)['b']*d$MedianAgeMarriage.s
# compute residual for each State
m.resid <- d$Marriage.s - mu

## R code 5.8
plot( Marriage.s ~ MedianAgeMarriage.s , d , col=rangi2 )
abline( m5.4 )
# loop over States
for ( i in 1:length(m.resid) ) {
    x <- d$MedianAgeMarriage.s[i] # x location of line segment
    y <- d$Marriage.s[i] # observed endpoint of line segment
    # draw the line segment
    lines( c(x,x) , c(mu[i],y) , lwd=0.5 , col=col.alpha("black",0.7) )
}

## R code 5.9
# prepare new counterfactual data
A.avg <- mean( d$MedianAgeMarriage.s )
R.seq <- seq( from=-3 , to=3 , length.out=30 )
pred.data <- data.frame(
    Marriage.s=R.seq,
    MedianAgeMarriage.s=A.avg
)

# compute counterfactual mean divorce (mu)
mu <- link( m5.3 , data=pred.data )
mu.mean <- apply( mu , 2 , mean )
mu.PI <- apply( mu , 2 , PI )

# simulate counterfactual divorce outcomes
R.sim <- sim( m5.3 , data=pred.data , n=1e4 )
R.PI <- apply( R.sim , 2 , PI )

# display predictions, hiding raw data with type="n"
plot( Divorce ~ Marriage.s , data=d , type="n" )
mtext( "MedianAgeMarriage.s = 0" )
lines( R.seq , mu.mean )
shade( mu.PI , R.seq )
shade( R.PI , R.seq )

## R code 5.10
R.avg <- mean( d$Marriage.s )
A.seq <- seq( from=-3 , to=3.5 , length.out=30 )
pred.data2 <- data.frame(
    Marriage.s=R.avg,
    MedianAgeMarriage.s=A.seq
)

mu <- link( m5.3 , data=pred.data2 )
mu.mean <- apply( mu , 2 , mean )
mu.PI <- apply( mu , 2 , PI )

A.sim <- sim( m5.3 , data=pred.data2 , n=1e4 )
A.PI <- apply( A.sim , 2 , PI )

plot( Divorce ~ MedianAgeMarriage.s , data=d , type="n" )
mtext( "Marriage.s = 0" )
lines( A.seq , mu.mean )
shade( mu.PI , A.seq )
shade( A.PI , A.seq )

## R code 5.11
# call link without specifying new data
# so it uses original data
mu <- link( m5.3 )

# summarize samples across cases
mu.mean <- apply( mu , 2 , mean )
mu.PI <- apply( mu , 2 , PI )

# simulate observations
# again no new data, so uses original data
divorce.sim <- sim( m5.3 , n=1e4 )
divorce.PI <- apply( divorce.sim , 2 , PI )

## R code 5.12
plot( mu.mean ~ d$Divorce , col=rangi2 , ylim=range(mu.PI) ,
    xlab="Observed divorce" , ylab="Predicted divorce" )
abline( a=0 , b=1 , lty=2 )
for ( i in 1:nrow(d) )
    lines( rep(d$Divorce[i],2) , c(mu.PI[1,i],mu.PI[2,i]) ,
        col=rangi2 )

## R code 5.13
identify( x=d$Divorce , y=mu.mean , labels=d$Loc , cex=0.8 )

## R code 5.14
# compute residuals
divorce.resid <- d$Divorce - mu.mean
# get ordering by divorce rate
o <- order(divorce.resid)
# make the plot
dotchart( divorce.resid[o] , labels=d$Loc[o] , xlim=c(-6,5) , cex=0.6 )
abline( v=0 , col=col.alpha("black",0.2) )
for ( i in 1:nrow(d) ) {
    j <- o[i] # which State in order
    lines( d$Divorce[j]-c(mu.PI[1,j],mu.PI[2,j]) , rep(i,2) )
    points( d$Divorce[j]-c(divorce.PI[1,j],divorce.PI[2,j]) , rep(i,2),
        pch=3 , cex=0.6 , col="gray" )
}

## R code 5.15
N <- 100                         # number of cases
x_real <- rnorm( N )             # x_real as Gaussian with mean 0 and stddev 1
x_spur <- rnorm( N , x_real )    # x_spur as Gaussian with mean=x_real
y <- rnorm( N , x_real )         # y as Gaussian with mean=x_real
d <- data.frame(y,x_real,x_spur) # bind all together in data frame

## R code 5.16
library(rethinking)
data(milk)
d <- milk
str(d)

## R code 5.17
m5.5 <- map(
    alist(
        kcal.per.g ~ dnorm( mu , sigma ) ,
        mu <- a + bn*neocortex.perc ,
        a ~ dnorm( 0 , 100 ) ,
        bn ~ dnorm( 0 , 1 ) ,
        sigma ~ dunif( 0 , 1 )
    ) ,
    data=d )

## R code 5.18
d$neocortex.perc

## R code 5.19
dcc <- d[ complete.cases(d) , ]

## R code 5.20
m5.5 <- map(
    alist(
        kcal.per.g ~ dnorm( mu , sigma ) ,
        mu <- a + bn*neocortex.perc ,
        a ~ dnorm( 0 , 100 ) ,
        bn ~ dnorm( 0 , 1 ) ,
        sigma ~ dunif( 0 , 1 )
    ) ,
    data=dcc )

## R code 5.21
precis( m5.5 , digits=3 )

## R code 5.22
coef(m5.5)["bn"] * ( 76 - 55 )

## R code 5.23
np.seq <- 0:100
pred.data <- data.frame( neocortex.perc=np.seq )

mu <- link( m5.5 , data=pred.data , n=1e4 )
mu.mean <- apply( mu , 2 , mean )
mu.PI <- apply( mu , 2 , PI )

plot( kcal.per.g ~ neocortex.perc , data=dcc , col=rangi2 )
lines( np.seq , mu.mean )
lines( np.seq , mu.PI[1,] , lty=2 )
lines( np.seq , mu.PI[2,] , lty=2 )

## R code 5.24
dcc$log.mass <- log(dcc$mass)

## R code 5.25
m5.6 <- map(
    alist(
        kcal.per.g ~ dnorm( mu , sigma ) ,
        mu <- a + bm*log.mass ,
        a ~ dnorm( 0 , 100 ) ,
        bm ~ dnorm( 0 , 1 ) ,
        sigma ~ dunif( 0 , 1 )
    ) ,
    data=dcc )
precis(m5.6)

## R code 5.26
m5.7 <- map(
    alist(
        kcal.per.g ~ dnorm( mu , sigma ) ,
        mu <- a + bn*neocortex.perc + bm*log.mass ,
        a ~ dnorm( 0 , 100 ) ,
        bn ~ dnorm( 0 , 1 ) ,
        bm ~ dnorm( 0 , 1 ) ,
        sigma ~ dunif( 0 , 1 )
    ) ,
    data=dcc )
precis(m5.7)

## R code 5.27
mean.log.mass <- mean( log(dcc$mass) )
np.seq <- 0:100
pred.data <- data.frame(
    neocortex.perc=np.seq,
    log.mass=mean.log.mass
)

mu <- link( m5.7 , data=pred.data , n=1e4 )
mu.mean <- apply( mu , 2 , mean )
mu.PI <- apply( mu , 2 , PI )

plot( kcal.per.g ~ neocortex.perc , data=dcc , type="n" )
lines( np.seq , mu.mean )
lines( np.seq , mu.PI[1,] , lty=2 )
lines( np.seq , mu.PI[2,] , lty=2 )

## R code 5.28
N <- 100                         # number of cases
rho <- 0.7                       # correlation btw x_pos and x_neg
x_pos <- rnorm( N )              # x_pos as Gaussian
x_neg <- rnorm( N , rho*x_pos ,  # x_neg correlated with x_pos
    sqrt(1-rho^2) )
y <- rnorm( N , x_pos - x_neg )  # y equally associated with x_pos, x_neg
d <- data.frame(y,x_pos,x_neg)   # bind all together in data frame

## R code 5.29
N <- 100                          # number of individuals
height <- rnorm(N,10,2)           # sim total height of each
leg_prop <- runif(N,0.4,0.5)      # leg as proportion of height
leg_left <- leg_prop*height +     # sim left leg as proportion + error
    rnorm( N , 0 , 0.02 )
leg_right <- leg_prop*height +    # sim right leg as proportion + error
    rnorm( N , 0 , 0.02 )
                                  # combine into data frame
d <- data.frame(height,leg_left,leg_right)

## R code 5.30
m5.8 <- map(
    alist(
        height ~ dnorm( mu , sigma ) ,
        mu <- a + bl*leg_left + br*leg_right ,
        a ~ dnorm( 10 , 100 ) ,
        bl ~ dnorm( 2 , 10 ) ,
        br ~ dnorm( 2 , 10 ) ,
        sigma ~ dunif( 0 , 10 )
    ) ,
    data=d )
precis(m5.8)

## R code 5.31
plot(precis(m5.8))

## R code 5.32
post <- extract.samples(m5.8)
plot( bl ~ br , post , col=col.alpha(rangi2,0.1) , pch=16 )

## R code 5.33
sum_blbr <- post$bl + post$br
dens( sum_blbr , col=rangi2 , lwd=2 , xlab="sum of bl and br" )

## R code 5.34
m5.9 <- map(
    alist(
        height ~ dnorm( mu , sigma ) ,
        mu <- a + bl*leg_left,
        a ~ dnorm( 10 , 100 ) ,
        bl ~ dnorm( 2 , 10 ) ,
        sigma ~ dunif( 0 , 10 )
    ) ,
    data=d )
precis(m5.9)

## R code 5.35
library(rethinking)
data(milk)
d <- milk

## R code 5.36
# kcal.per.g regressed on perc.fat
m5.10 <- map(
    alist(
        kcal.per.g ~ dnorm( mu , sigma ) ,
        mu <- a + bf*perc.fat ,
        a ~ dnorm( 0.6 , 10 ) ,
        bf ~ dnorm( 0 , 1 ) ,
        sigma ~ dunif( 0 , 10 )
    ) ,
    data=d )

# kcal.per.g regressed on perc.lactose
m5.11 <- map(
    alist(
        kcal.per.g ~ dnorm( mu , sigma ) ,
        mu <- a + bl*perc.lactose ,
        a ~ dnorm( 0.6 , 10 ) ,
        bl ~ dnorm( 0 , 1 ) ,
        sigma ~ dunif( 0 , 10 )
    ) ,
    data=d )

precis( m5.10 , digits=3 )
precis( m5.11 , digits=3 )

## R code 5.37
m5.12 <- map(
    alist(
        kcal.per.g ~ dnorm( mu , sigma ) ,
        mu <- a + bf*perc.fat + bl*perc.lactose ,
        a ~ dnorm( 0.6 , 10 ) ,
        bf ~ dnorm( 0 , 1 ) ,
        bl ~ dnorm( 0 , 1 ) ,
        sigma ~ dunif( 0 , 10 )
    ) ,
    data=d )
precis( m5.12 , digits=3 )

## R code 5.38
pairs( ~ kcal.per.g + perc.fat + perc.lactose ,
    data=d , col=rangi2 )

## R code 5.39
cor( d$perc.fat , d$perc.lactose )

## R code 5.40
library(rethinking)
data(milk)
d <- milk
sim.coll <- function( r=0.9 ) {
    d$x <- rnorm( nrow(d) , mean=r*d$perc.fat ,
        sd=sqrt( (1-r^2)*var(d$perc.fat) ) )
    m <- lm( kcal.per.g ~ perc.fat + x , data=d )
    sqrt( diag( vcov(m) ) )[2] # stddev of parameter
}
rep.sim.coll <- function( r=0.9 , n=100 ) {
    stddev <- replicate( n , sim.coll(r) )
    mean(stddev)
}
r.seq <- seq(from=0,to=0.99,by=0.01)
stddev <- sapply( r.seq , function(z) rep.sim.coll(r=z,n=100) )
plot( stddev ~ r.seq , type="l" , col=rangi2, lwd=2 , xlab="correlation" )

## R code 5.41
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

## R code 5.42
m5.13 <- map(
    alist(
        h1 ~ dnorm(mu,sigma),
        mu <- a + bh*h0 + bt*treatment + bf*fungus,
        a ~ dnorm(0,100),
        c(bh,bt,bf) ~ dnorm(0,10),
        sigma ~ dunif(0,10)
    ),
    data=d )
precis(m5.13)

## R code 5.43
m5.14 <- map(
    alist(
        h1 ~ dnorm(mu,sigma),
        mu <- a + bh*h0 + bt*treatment,
        a ~ dnorm(0,100),
        c(bh,bt) ~ dnorm(0,10),
        sigma ~ dunif(0,10)
    ),
    data=d )
precis(m5.14)

## R code 5.44
data(Howell1)
d <- Howell1
str(d)

## R code 5.45
m5.15 <- map(
    alist(
        height ~ dnorm( mu , sigma ) ,
        mu <- a + bm*male ,
        a ~ dnorm( 178 , 100 ) ,
        bm ~ dnorm( 0 , 10 ) ,
        sigma ~ dunif( 0 , 50 )
    ) ,
    data=d )
precis(m5.15)

## R code 5.46
post <- extract.samples(m5.15)
mu.male <- post$a + post$bm
PI(mu.male)

## R code 5.47
m5.15b <- map(
    alist(
        height ~ dnorm( mu , sigma ) ,
        mu <- af*(1-male) + am*male ,
        af ~ dnorm( 178 , 100 ) ,
        am ~ dnorm( 178 , 100 ) ,
        sigma ~ dunif( 0 , 50 )
    ) ,
    data=d )

## R code 5.48
data(milk)
d <- milk
unique(d$clade)

## R code 5.49
( d$clade.NWM <- ifelse( d$clade=="New World Monkey" , 1 , 0 ) )

## R code 5.50
d$clade.OWM <- ifelse( d$clade=="Old World Monkey" , 1 , 0 )
d$clade.S <- ifelse( d$clade=="Strepsirrhine" , 1 , 0 )

## R code 5.51
m5.16 <- map(
    alist(
        kcal.per.g ~ dnorm( mu , sigma ) ,
        mu <- a + b.NWM*clade.NWM + b.OWM*clade.OWM + b.S*clade.S ,
        a ~ dnorm( 0.6 , 10 ) ,
        b.NWM ~ dnorm( 0 , 1 ) ,
        b.OWM ~ dnorm( 0 , 1 ) ,
        b.S ~ dnorm( 0 , 1 ) ,
        sigma ~ dunif( 0 , 10 )
    ) ,
    data=d )
precis(m5.16)

## R code 5.52
# sample posterior
post <- extract.samples(m5.16)

# compute averages for each category
mu.ape <- post$a
mu.NWM <- post$a + post$b.NWM
mu.OWM <- post$a + post$b.OWM
mu.S <- post$a + post$b.S

# summarize using precis
precis( data.frame(mu.ape,mu.NWM,mu.OWM,mu.S) )

## R code 5.53
diff.NWM.OWM <- mu.NWM - mu.OWM
quantile( diff.NWM.OWM , probs=c(0.025,0.5,0.975) )

## R code 5.54
( d$clade_id <- coerce_index(d$clade) )

## R code 5.55
m5.16_alt <- map(
    alist(
        kcal.per.g ~ dnorm( mu , sigma ) ,
        mu <- a[clade_id] ,
        a[clade_id] ~ dnorm( 0.6 , 10 ) ,
        sigma ~ dunif( 0 , 10 )
    ) ,
    data=d )
precis( m5.16_alt , depth=2 )

## R code 5.56
m5.17 <- lm( y ~ 1 + x , data=d )
m5.18 <- lm( y ~ 1 + x + z + w , data=d )

## R code 5.57
m5.17 <- lm( y ~ 1 + x , data=d )
m5.19 <- lm( y ~ x , data=d )

## R code 5.58
m5.20 <- lm( y ~ 0 + x , data=d )
m5.21 <- lm( y ~ x - 1 , data=d )

## R code 5.59
m5.22 <- lm( y ~ 1 + as.factor(season) , data=d )

## R code 5.60
d$x2 <- d$x^2
d$x3 <- d$x^3
m5.23 <- lm( y ~ 1 + x + x2 + x3 , data=d )

## R code 5.61
m5.24 <- lm( y ~ 1 + x + I(x^2) + I(x^3) , data=d )

## R code 5.62
data(cars)
glimmer( dist ~ speed , data=cars )

## R code 6.1
sppnames <- c( "afarensis","africanus","habilis","boisei",
    "rudolfensis","ergaster","sapiens")
brainvolcc <- c( 438 , 452 , 612, 521, 752, 871, 1350 )
masskg <- c( 37.0 , 35.5 , 34.5 , 41.5 , 55.5 , 61.0 , 53.5 )
d <- data.frame( species=sppnames , brain=brainvolcc , mass=masskg )

## R code 6.2
m6.1 <- lm( brain ~ mass , data=d )

## R code 6.3
1 - var(resid(m6.1))/var(d$brain)

## R code 6.4
m6.2 <- lm( brain ~ mass + I(mass^2) , data=d )

## R code 6.5
m6.3 <- lm( brain ~ mass + I(mass^2) + I(mass^3) , data=d )
m6.4 <- lm( brain ~ mass + I(mass^2) + I(mass^3) + I(mass^4) ,
    data=d )
m6.5 <- lm( brain ~ mass + I(mass^2) + I(mass^3) + I(mass^4) +
    I(mass^5) , data=d )
m6.6 <- lm( brain ~ mass + I(mass^2) + I(mass^3) + I(mass^4) +
    I(mass^5) + I(mass^6) , data=d )

## R code 6.6
m6.7 <- lm( brain ~ 1 , data=d )

## R code 6.7
d.new <- d[ -i , ]

## R code 6.8
plot( brain ~ mass , d , col="slateblue" )
for ( i in 1:nrow(d) ) {
    d.new <- d[ -i , ]
    m0 <- lm( brain ~ mass, d.new )
    abline( m0 , col=col.alpha("black",0.5) )
}

## R code 6.9
p <- c( 0.3 , 0.7 )
-sum( p*log(p) )

## R code 6.10
# fit model with lm
m6.1 <- lm( brain ~ mass , d )

# compute deviance by cheating
(-2) * logLik(m6.1)

## R code 6.11
# standardize the mass before fitting
d$mass.s <- (d$mass-mean(d$mass))/sd(d$mass)
m6.8 <- map(
    alist(
        brain ~ dnorm( mu , sigma ) ,
        mu <- a + b*mass.s
    ) ,
    data=d ,
    start=list(a=mean(d$brain),b=0,sigma=sd(d$brain)) ,
    method="Nelder-Mead" )

# extract MAP estimates
theta <- coef(m6.8)

# compute deviance
dev <- (-2)*sum( dnorm(
            d$brain ,
            mean=theta[1]+theta[2]*d$mass.s ,
            sd=theta[3] ,
            log=TRUE ) )
dev

## R code 6.12
N <- 20
kseq <- 1:5
dev <- sapply( kseq , function(k) {
        print(k);
        r <- replicate( 1e4 , sim.train.test( N=N, k=k ) );
        c( mean(r[1,]) , mean(r[2,]) , sd(r[1,]) , sd(r[2,]) )
    } )

## R code 6.13
        r <- mcreplicate( 1e4 , sim.train.test( N=N, k=k ) , mc.cores=4 )

## R code 6.14
plot( 1:5 , dev[1,] , ylim=c( min(dev[1:2,])-5 , max(dev[1:2,])+10 ) ,
    xlim=c(1,5.1) , xlab="number of parameters" , ylab="deviance" ,
    pch=16 , col=rangi2 )
mtext( concat( "N = ",N ) )
points( (1:5)+0.1 , dev[2,] )
for ( i in kseq ) {
    pts_in <- dev[1,i] + c(-1,+1)*dev[3,i]
    pts_out <- dev[2,i] + c(-1,+1)*dev[4,i]
    lines( c(i,i) , pts_in , col=rangi2 )
    lines( c(i,i)+0.1 , pts_out )
}

## R code 6.15
data(cars)
m <- map(
    alist(
        dist ~ dnorm(mu,sigma),
        mu <- a + b*speed,
        a ~ dnorm(0,100),
        b ~ dnorm(0,10),
        sigma ~ dunif(0,30)
    ) , data=cars )
post <- extract.samples(m,n=1000)

## R code 6.16
n_samples <- 1000
ll <- sapply( 1:n_samples ,
    function(s) {
        mu <- post$a[s] + post$b[s]*cars$speed
        dnorm( cars$dist , mu , post$sigma[s] , log=TRUE )
    } )

## R code 6.17
n_cases <- nrow(cars)
lppd <- sapply( 1:n_cases , function(i) log_sum_exp(ll[i,]) - log(n_samples) )

## R code 6.18
pWAIC <- sapply( 1:n_cases , function(i) var(ll[i,]) )

## R code 6.19
-2*( sum(lppd) - sum(pWAIC) )

## R code 6.20
waic_vec <- -2*( lppd - pWAIC )
sqrt( n_cases*var(waic_vec) )

## R code 6.21
data(milk)
d <- milk[ complete.cases(milk) , ]
d$neocortex <- d$neocortex.perc / 100
dim(d)

## R code 6.22
a.start <- mean(d$kcal.per.g)
sigma.start <- log(sd(d$kcal.per.g))
m6.11 <- map(
    alist(
        kcal.per.g ~ dnorm( a , exp(log.sigma) )
    ) ,
    data=d , start=list(a=a.start,log.sigma=sigma.start) )
m6.12 <- map(
    alist(
        kcal.per.g ~ dnorm( mu , exp(log.sigma) ) ,
        mu <- a + bn*neocortex
    ) ,
    data=d , start=list(a=a.start,bn=0,log.sigma=sigma.start) )
m6.13 <- map(
    alist(
        kcal.per.g ~ dnorm( mu , exp(log.sigma) ) ,
        mu <- a + bm*log(mass)
    ) ,
    data=d , start=list(a=a.start,bm=0,log.sigma=sigma.start) )
m6.14 <- map(
    alist(
        kcal.per.g ~ dnorm( mu , exp(log.sigma) ) ,
        mu <- a + bn*neocortex + bm*log(mass)
    ) ,
    data=d , start=list(a=a.start,bn=0,bm=0,log.sigma=sigma.start) )

## R code 6.23
WAIC( m6.14 )

## R code 6.24
( milk.models <- compare( m6.11 , m6.12 , m6.13 , m6.14 ) )

## R code 6.25
plot( milk.models , SE=TRUE , dSE=TRUE )

## R code 6.26
diff <- rnorm( 1e5 , 6.7 , 7.26 )
sum(diff<0)/1e5

## R code 6.27
coeftab(m6.11,m6.12,m6.13,m6.14)

## R code 6.28
plot( coeftab(m6.11,m6.12,m6.13,m6.14) )

## R code 6.29
# compute counterfactual predictions
# neocortex from 0.5 to 0.8
nc.seq <- seq(from=0.5,to=0.8,length.out=30)
d.predict <- list(
    kcal.per.g = rep(0,30), # empty outcome
    neocortex = nc.seq,     # sequence of neocortex
    mass = rep(4.5,30)      # average mass
)
pred.m6.14 <- link( m6.14 , data=d.predict )
mu <- apply( pred.m6.14 , 2 , mean )
mu.PI <- apply( pred.m6.14 , 2 , PI )

# plot it all
plot( kcal.per.g ~ neocortex , d , col=rangi2 )
lines( nc.seq , mu , lty=2 )
lines( nc.seq , mu.PI[1,] , lty=2 )
lines( nc.seq , mu.PI[2,] , lty=2 )

## R code 6.30
milk.ensemble <- ensemble( m6.11 , m6.12 , m6.13 , m6.14 , data=d.predict )
mu <- apply( milk.ensemble$link , 2 , mean )
mu.PI <- apply( milk.ensemble$link , 2 , PI )
lines( nc.seq , mu )
shade( mu.PI , nc.seq )

## R code 6.31
library(rethinking)
data(Howell1)
d <- Howell1
d$age <- (d$age - mean(d$age))/sd(d$age)
set.seed( 1000 )
i <- sample(1:nrow(d),size=nrow(d)/2)
d1 <- d[ i , ]
d2 <- d[ -i , ]

## R code 6.32
sum( dnorm( d2$height , mu , sigma , log=TRUE ) )

## R code 7.1
library(rethinking)
data(rugged)
d <- rugged

# make log version of outcome
d$log_gdp <- log( d$rgdppc_2000 )

# extract countries with GDP data
dd <- d[ complete.cases(d$rgdppc_2000) , ]

# split countries into Africa and not-Africa
d.A1 <- dd[ dd$cont_africa==1 , ] # Africa
d.A0 <- dd[ dd$cont_africa==0 , ] # not Africa

## R code 7.2
# African nations
m7.1 <- map(
    alist(
        log_gdp ~ dnorm( mu , sigma ) ,
        mu <- a + bR*rugged ,
        a ~ dnorm( 8 , 100 ) ,
        bR ~ dnorm( 0 , 1 ) ,
        sigma ~ dunif( 0 , 10 )
    ) ,
    data=d.A1 )

# non-African nations
m7.2 <- map(
    alist(
        log_gdp ~ dnorm( mu , sigma ) ,
        mu <- a + bR*rugged ,
        a ~ dnorm( 8 , 100 ) ,
        bR ~ dnorm( 0 , 1 ) ,
        sigma ~ dunif( 0 , 10 )
    ) ,
    data=d.A0 )

## R code 7.3
m7.3 <- map(
    alist(
        log_gdp ~ dnorm( mu , sigma ) ,
        mu <- a + bR*rugged ,
        a ~ dnorm( 8 , 100 ) ,
        bR ~ dnorm( 0 , 1 ) ,
        sigma ~ dunif( 0 , 10 )
    ) ,
    data=dd )

## R code 7.4
m7.4 <- map(
    alist(
        log_gdp ~ dnorm( mu , sigma ) ,
        mu <- a + bR*rugged + bA*cont_africa ,
        a ~ dnorm( 8 , 100 ) ,
        bR ~ dnorm( 0 , 1 ) ,
        bA ~ dnorm( 0 , 1 ) ,
        sigma ~ dunif( 0 , 10 )
    ) ,
    data=dd )

## R code 7.5
compare( m7.3 , m7.4 )

## R code 7.6
rugged.seq <- seq(from=-1,to=8,by=0.25)

# compute mu over samples, fixing cont_africa=0
mu.NotAfrica <- link( m7.4 , data=data.frame(cont_africa=0,rugged=rugged.seq) )

# compute mu over samples, fixing cont_africa=1
mu.Africa <- link( m7.4 , data=data.frame(cont_africa=1,rugged=rugged.seq) )

# summarize to means and intervals
mu.NotAfrica.mean <- apply( mu.NotAfrica , 2 , mean )
mu.NotAfrica.PI <- apply( mu.NotAfrica , 2 , PI , prob=0.97 )
mu.Africa.mean <- apply( mu.Africa , 2 , mean )
mu.Africa.PI <- apply( mu.Africa , 2 , PI , prob=0.97 )

## R code 7.7
m7.5 <- map(
    alist(
        log_gdp ~ dnorm( mu , sigma ) ,
        mu <- a + gamma*rugged + bA*cont_africa ,
        gamma <- bR + bAR*cont_africa ,
        a ~ dnorm( 8 , 100 ) ,
        bA ~ dnorm( 0 , 1 ) ,
        bR ~ dnorm( 0 , 1 ) ,
        bAR ~ dnorm( 0 , 1 ) ,
        sigma ~ dunif( 0 , 10 )
    ) ,
    data=dd )

## R code 7.8
compare( m7.3 , m7.4 , m7.5 )

## R code 7.9
m7.5b <- map(
    alist(
        log_gdp ~ dnorm( mu , sigma ) ,
        mu <- a + bR*rugged + bAR*rugged*cont_africa + bA*cont_africa,
        a ~ dnorm( 8 , 100 ) ,
        bA ~ dnorm( 0 , 1 ) ,
        bR ~ dnorm( 0 , 1 ) ,
        bAR ~ dnorm( 0 , 1 ) ,
        sigma ~ dunif( 0 , 10 )
    ) ,
    data=dd )

## R code 7.10
rugged.seq <- seq(from=-1,to=8,by=0.25)

mu.Africa <- link( m7.5 , data=data.frame(cont_africa=1,rugged=rugged.seq) )
mu.Africa.mean <- apply( mu.Africa , 2 , mean )
mu.Africa.PI <- apply( mu.Africa , 2 , PI , prob=0.97 )

mu.NotAfrica <- link( m7.5 , data=data.frame(cont_africa=0,rugged=rugged.seq) )
mu.NotAfrica.mean <- apply( mu.NotAfrica , 2 , mean )
mu.NotAfrica.PI <- apply( mu.NotAfrica , 2 , PI , prob=0.97 )

## R code 7.11
# plot African nations with regression
d.A1 <- dd[dd$cont_africa==1,]
plot( log(rgdppc_2000) ~ rugged , data=d.A1 ,
    col=rangi2 , ylab="log GDP year 2000" ,
    xlab="Terrain Ruggedness Index" )
mtext( "African nations" , 3 )
lines( rugged.seq , mu.Africa.mean , col=rangi2 )
shade( mu.Africa.PI , rugged.seq , col=col.alpha(rangi2,0.3) )

# plot non-African nations with regression
d.A0 <- dd[dd$cont_africa==0,]
plot( log(rgdppc_2000) ~ rugged , data=d.A0 ,
    col="black" , ylab="log GDP year 2000" ,
    xlab="Terrain Ruggedness Index" )
mtext( "Non-African nations" , 3 )
lines( rugged.seq , mu.NotAfrica.mean )
shade( mu.NotAfrica.PI , rugged.seq )

## R code 7.12
precis(m7.5)

## R code 7.13
post <- extract.samples( m7.5 )
gamma.Africa <- post$bR + post$bAR*1
gamma.notAfrica <- post$bR + post$bAR*0

## R code 7.14
mean( gamma.Africa)
mean( gamma.notAfrica )

## R code 7.15
dens( gamma.Africa , xlim=c(-0.5,0.6) , ylim=c(0,5.5) ,
    xlab="gamma" , col=rangi2 )
dens( gamma.notAfrica , add=TRUE )

## R code 7.16
diff <- gamma.Africa - gamma.notAfrica
sum( diff < 0 ) / length( diff )

## R code 7.17
# get minimum and maximum rugged values
q.rugged <- range(dd$rugged)

# compute lines and confidence intervals
mu.ruggedlo <- link( m7.5 ,
    data=data.frame(rugged=q.rugged[1],cont_africa=0:1) )
mu.ruggedlo.mean <- apply( mu.ruggedlo , 2 , mean )
mu.ruggedlo.PI <- apply( mu.ruggedlo , 2 , PI )

mu.ruggedhi <- link( m7.5 ,
    data=data.frame(rugged=q.rugged[2],cont_africa=0:1) )
mu.ruggedhi.mean <- apply( mu.ruggedhi , 2 , mean )
mu.ruggedhi.PI <- apply( mu.ruggedhi , 2 , PI )

# plot it all, splitting points at median
med.r <- median(dd$rugged)
ox <- ifelse( dd$rugged > med.r , 0.05 , -0.05 )
plot( dd$cont_africa + ox , log(dd$rgdppc_2000) ,
    col=ifelse(dd$rugged>med.r,rangi2,"black") ,
    xlim=c(-0.25,1.25) , xaxt="n" , ylab="log GDP year 2000" ,
    xlab="Continent" )
axis( 1 , at=c(0,1) , labels=c("other","Africa") )
lines( 0:1 , mu.ruggedlo.mean , lty=2 )
shade( mu.ruggedlo.PI , 0:1 )
lines( 0:1 , mu.ruggedhi.mean , col=rangi2 )
shade( mu.ruggedhi.PI , 0:1 , col=col.alpha(rangi2,0.25) )

## R code 7.18
library(rethinking)
data(tulips)
d <- tulips
str(d)

## R code 7.19
m7.6 <- map(
    alist(
        blooms ~ dnorm( mu , sigma ) ,
        mu <- a + bW*water + bS*shade ,
        a ~ dnorm( 0 , 100 ) ,
        bW ~ dnorm( 0 , 100 ) ,
        bS ~ dnorm( 0 , 100 ) ,
        sigma ~ dunif( 0 , 100 )
    ) ,
    data=d )
m7.7 <- map(
    alist(
        blooms ~ dnorm( mu , sigma ) ,
        mu <- a + bW*water + bS*shade + bWS*water*shade ,
        a ~ dnorm( 0 , 100 ) ,
        bW ~ dnorm( 0 , 100 ) ,
        bS ~ dnorm( 0 , 100 ) ,
        bWS ~ dnorm( 0 , 100 ) ,
        sigma ~ dunif( 0 , 100 )
    ) ,
    data=d )

## R code 7.20
m7.6 <- map(
    alist(
        blooms ~ dnorm( mu , sigma ) ,
        mu <- a + bW*water + bS*shade ,
        a ~ dnorm( 0 , 100 ) ,
        bW ~ dnorm( 0 , 100 ) ,
        bS ~ dnorm( 0 , 100 ) ,
        sigma ~ dunif( 0 , 100 )
    ) ,
    data=d ,
    method="Nelder-Mead" ,
    control=list(maxit=1e4) )
m7.7 <- map(
    alist(
        blooms ~ dnorm( mu , sigma ) ,
        mu <- a + bW*water + bS*shade + bWS*water*shade ,
        a ~ dnorm( 0 , 100 ) ,
        bW ~ dnorm( 0 , 100 ) ,
        bS ~ dnorm( 0 , 100 ) ,
        bWS ~ dnorm( 0 , 100 ) ,
        sigma ~ dunif( 0 , 100 )
    ) ,
    data=d ,
    method="Nelder-Mead" ,
    control=list(maxit=1e4) )

## R code 7.21
coeftab(m7.6,m7.7)

## R code 7.22
compare( m7.6 , m7.7 )

## R code 7.23
d$shade.c <- d$shade - mean(d$shade)
d$water.c <- d$water - mean(d$water)

## R code 7.24
m7.8 <- map(
    alist(
        blooms ~ dnorm( mu , sigma ) ,
        mu <- a + bW*water.c + bS*shade.c ,
        a ~ dnorm( 130 , 100 ) ,
        bW ~ dnorm( 0 , 100 ) ,
        bS ~ dnorm( 0 , 100 ) ,
        sigma ~ dunif( 0 , 100 )
    ) ,
    data=d ,
    start=list(a=mean(d$blooms),bW=0,bS=0,sigma=sd(d$blooms)) )
m7.9 <- map(
    alist(
        blooms ~ dnorm( mu , sigma ) ,
        mu <- a + bW*water.c + bS*shade.c + bWS*water.c*shade.c ,
        a ~ dnorm( 130 , 100 ) ,
        bW ~ dnorm( 0 , 100 ) ,
        bS ~ dnorm( 0 , 100 ) ,
        bWS ~ dnorm( 0 , 100 ) ,
        sigma ~ dunif( 0 , 100 )
    ) ,
    data=d ,
    start=list(a=mean(d$blooms),bW=0,bS=0,bWS=0,sigma=sd(d$blooms)) )
coeftab(m7.8,m7.9)

## R code 7.25
k <- coef(m7.7)
k[1] + k[2]*2 + k[3]*2 + k[4]*2*2

## R code 7.26
k <- coef(m7.9)
k[1] + k[2]*0 + k[3]*0 + k[4]*0*0

## R code 7.27
precis(m7.9)

## R code 7.28
# make a plot window with three panels in a single row
par(mfrow=c(1,3)) # 1 row, 3 columns

# loop over values of water.c and plot predictions
shade.seq <- -1:1
for ( w in -1:1 ) {
    dt <- d[d$water.c==w,]
    plot( blooms ~ shade.c , data=dt , col=rangi2 ,
        main=paste("water.c =",w) , xaxp=c(-1,1,2) , ylim=c(0,362) ,
        xlab="shade (centered)" )
    mu <- link( m7.9 , data=data.frame(water.c=w,shade.c=shade.seq) )
    mu.mean <- apply( mu , 2 , mean )
    mu.PI <- apply( mu , 2 , PI , prob=0.97 )
    lines( shade.seq , mu.mean )
    lines( shade.seq , mu.PI[1,] , lty=2 )
    lines( shade.seq , mu.PI[2,] , lty=2 )
}

## R code 7.29
m7.x <- lm( y ~ x + z + x*z , data=d )

## R code 7.30
m7.x <- lm( y ~ x*z , data=d )

## R code 7.31
m7.x <- lm( y ~ x + x*z - z , data=d )

## R code 7.32
m7.x <- lm( y ~ x*z*w , data=d )

## R code 7.33
x <- z <- w <- 1
colnames( model.matrix(~x*z*w) )

## R code 7.34
d$lang.per.cap <- d$num.lang / d$k.pop

## R code 8.1
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

## R code 8.2
library(rethinking)
data(rugged)
d <- rugged
d$log_gdp <- log(d$rgdppc_2000)
dd <- d[ complete.cases(d$rgdppc_2000) , ]

## R code 8.3
m8.1 <- map(
    alist(
        log_gdp ~ dnorm( mu , sigma ) ,
        mu <- a + bR*rugged + bA*cont_africa + bAR*rugged*cont_africa ,
        a ~ dnorm(0,100),
        bR ~ dnorm(0,10),
        bA ~ dnorm(0,10),
        bAR ~ dnorm(0,10),
        sigma ~ dunif(0,10)
    ) ,
    data=dd )
precis(m8.1)

## R code 8.4
dd.trim <- dd[ , c("log_gdp","rugged","cont_africa") ]
str(dd.trim)

## R code 8.5
m8.1stan <- map2stan(
    alist(
        log_gdp ~ dnorm( mu , sigma ) ,
        mu <- a + bR*rugged + bA*cont_africa + bAR*rugged*cont_africa ,
        a ~ dnorm(0,100),
        bR ~ dnorm(0,10),
        bA ~ dnorm(0,10),
        bAR ~ dnorm(0,10),
        sigma ~ dcauchy(0,2)
    ) ,
    data=dd.trim )

## R code 8.6
precis(m8.1stan)

## R code 8.7
m8.1stan_4chains <- map2stan( m8.1stan , chains=4 , cores=4 )
precis(m8.1stan_4chains)

## R code 8.8
post <- extract.samples( m8.1stan )
str(post)

## R code 8.9
pairs(post)

## R code 8.10
pairs(m8.1stan)

## R code 8.11
show(m8.1stan)

## R code 8.12
plot(m8.1stan)

## R code 8.13
y <- c(-1,1)
m8.2 <- map2stan(
    alist(
        y ~ dnorm( mu , sigma ) ,
        mu <- alpha
    ) ,
    data=list(y=y) , start=list(alpha=0,sigma=1) ,
    chains=2 , iter=4000 , warmup=1000 )

## R code 8.14
precis(m8.2)

## R code 8.15
m8.3 <- map2stan(
    alist(
        y ~ dnorm( mu , sigma ) ,
        mu <- alpha ,
        alpha ~ dnorm( 1 , 10 ) ,
        sigma ~ dcauchy( 0 , 1 )
    ) ,
    data=list(y=y) , start=list(alpha=0,sigma=1) ,
    chains=2 , iter=4000 , warmup=1000 )
precis(m8.3)

## R code 8.16
y <- rcauchy(1e4,0,5)
mu <- sapply( 1:length(y) , function(i) sum(y[1:i])/i )
plot(mu,type="l")

## R code 8.17
y <- rnorm( 100 , mean=0 , sd=1 )

## R code 8.18
m8.4 <- map2stan(
    alist(
        y ~ dnorm( mu , sigma ) ,
        mu <- a1 + a2 ,
        sigma ~ dcauchy( 0 , 1 )
    ) ,
    data=list(y=y) , start=list(a1=0,a2=0,sigma=1) ,
    chains=2 , iter=4000 , warmup=1000 )
precis(m8.4)

## R code 8.19
m8.5 <- map2stan(
    alist(
        y ~ dnorm( mu , sigma ) ,
        mu <- a1 + a2 ,
        a1 ~ dnorm( 0 , 10 ) ,
        a2 ~ dnorm( 0 , 10 ) ,
        sigma ~ dcauchy( 0 , 1 )
    ) ,
    data=list(y=y) , start=list(a1=0,a2=0,sigma=1) ,
    chains=2 , iter=4000 , warmup=1000 )
precis(m8.5)

## R code 8.20
mp <- map2stan(
    alist(
        a ~ dnorm(0,1),
        b ~ dcauchy(0,1)
    ),
    data=list(y=1),
    start=list(a=0,b=0),
    iter=1e4, warmup=100 , WAIC=FALSE )

## R code 8.21
N <- 100                          # number of individuals
height <- rnorm(N,10,2)           # sim total height of each
leg_prop <- runif(N,0.4,0.5)      # leg as proportion of height
leg_left <- leg_prop*height +     # sim left leg as proportion + error
    rnorm( N , 0 , 0.02 )
leg_right <- leg_prop*height +    # sim right leg as proportion + error
    rnorm( N , 0 , 0.02 )
                                  # combine into data frame
d <- data.frame(height,leg_left,leg_right)

## R code 8.22
m5.8s <- map2stan(
    alist(
        height ~ dnorm( mu , sigma ) ,
        mu <- a + bl*leg_left + br*leg_right ,
        a ~ dnorm( 10 , 100 ) ,
        bl ~ dnorm( 2 , 10 ) ,
        br ~ dnorm( 2 , 10 ) ,
        sigma ~ dcauchy( 0 , 1 )
    ) ,
    data=d, chains=4,
    start=list(a=10,bl=0,br=0,sigma=1) )

## R code 8.23
m5.8s2 <- map2stan(
    alist(
        height ~ dnorm( mu , sigma ) ,
        mu <- a + bl*leg_left + br*leg_right ,
        a ~ dnorm( 10 , 100 ) ,
        bl ~ dnorm( 2 , 10 ) ,
        br ~ dnorm( 2 , 10 ) & T[0,] ,
        sigma ~ dcauchy( 0 , 1 )
    ) ,
    data=d, chains=4,
    start=list(a=10,bl=0,br=0,sigma=1) )

## R code 9.1
p <- list()
p$A <- c(0,0,10,0,0)
p$B <- c(0,1,8,1,0)
p$C <- c(0,2,6,2,0)
p$D <- c(1,2,4,2,1)
p$E <- c(2,2,2,2,2)

## R code 9.2
p_norm <- lapply( p , function(q) q/sum(q))

## R code 9.3
( H <- sapply( p_norm , function(q) -sum(ifelse(q==0,0,q*log(q))) ) )

## R code 9.4
ways <- c(1,90,1260,37800,113400)
logwayspp <- log(ways)/10

## R code 9.5
# build list of the candidate distributions
p <- list()
p[[1]] <- c(1/4,1/4,1/4,1/4)
p[[2]] <- c(2/6,1/6,1/6,2/6)
p[[3]] <- c(1/6,2/6,2/6,1/6)
p[[4]] <- c(1/8,4/8,2/8,1/8)

# compute expected value of each
sapply( p , function(p) sum(p*c(0,1,1,2)) )

## R code 9.6
# compute entropy of each distribution
sapply( p , function(p) -sum( p*log(p) ) )

## R code 9.7
p <- 0.7
( A <- c( (1-p)^2 , p*(1-p) , (1-p)*p , p^2 ) )

## R code 9.8
-sum( A*log(A) )

## R code 9.9
sim.p <- function(G=1.4) {
    x123 <- runif(3)
    x4 <- ( (G)*sum(x123)-x123[2]-x123[3] )/(2-G)
    z <- sum( c(x123,x4) )
    p <- c( x123 , x4 )/z
    list( H=-sum( p*log(p) ) , p=p )
}

## R code 9.10
H <- replicate( 1e5 , sim.p(1.4) )
dens( as.numeric(H[1,]) , adj=0.1 )

## R code 9.11
entropies <- as.numeric(H[1,])
distributions <- H[2,]

## R code 9.12
max(entropies)

## R code 9.13
distributions[ which.max(entropies) ]

## R code 10.1
library(rethinking)
data(chimpanzees)
d <- chimpanzees

## R code 10.2
m10.1 <- map(
    alist(
        pulled_left ~ dbinom( 1 , p ) ,
        logit(p) <- a ,
        a ~ dnorm(0,10)
    ) ,
    data=d )
precis(m10.1)

## R code 10.3
logistic( c(0.18,0.46) )

## R code 10.4
m10.2 <- map(
    alist(
        pulled_left ~ dbinom( 1 , p ) ,
        logit(p) <- a + bp*prosoc_left ,
        a ~ dnorm(0,10) ,
        bp ~ dnorm(0,10)
    ) ,
    data=d )
m10.3 <- map(
    alist(
        pulled_left ~ dbinom( 1 , p ) ,
        logit(p) <- a + (bp + bpC*condition)*prosoc_left ,
        a ~ dnorm(0,10) ,
        bp ~ dnorm(0,10) ,
        bpC ~ dnorm(0,10)
    ) ,
    data=d )

## R code 10.5
compare( m10.1 , m10.2 , m10.3 )

## R code 10.6
precis(m10.3)

## R code 10.7
exp(0.61)

## R code 10.8
logistic( 4 )

## R code 10.9
logistic( 4 + 0.61 )

## R code 10.10
# dummy data for predictions across treatments
d.pred <- data.frame(
    prosoc_left = c(0,1,0,1),   # right/left/right/left
    condition = c(0,0,1,1)      # control/control/partner/partner
)

# build prediction ensemble
chimp.ensemble <- ensemble( m10.1 , m10.2 , m10.3 , data=d.pred )

# summarize
pred.p <- apply( chimp.ensemble$link , 2 , mean )
pred.p.PI <- apply( chimp.ensemble$link , 2 , PI )

## R code 10.11
# empty plot frame with good axes
plot( 0 , 0 , type="n" , xlab="prosoc_left/condition" ,
    ylab="proportion pulled left" , ylim=c(0,1) , xaxt="n" ,
    xlim=c(1,4) )
axis( 1 , at=1:4 , labels=c("0/0","1/0","0/1","1/1") )

# plot raw data, one trend for each of 7 individual chimpanzees
# will use by() here; see Overthinking box for explanation
p <- by( d$pulled_left ,
    list(d$prosoc_left,d$condition,d$actor) , mean )
for ( chimp in 1:7 )
    lines( 1:4 , as.vector(p[,,chimp]) , col=rangi2 , lwd=1.5 )

# now superimpose posterior predictions
lines( 1:4 , pred.p )
shade( pred.p.PI , 1:4 )

## R code 10.12
# clean NAs from the data
d2 <- d
d2$recipient <- NULL

# re-use map fit to get the formula
m10.3stan <- map2stan( m10.3 , data=d2 , iter=1e4 , warmup=1000 )
precis(m10.3stan)

## R code 10.13
pairs(m10.3stan)

## R code 10.14
m10.4 <- map2stan(
    alist(
        pulled_left ~ dbinom( 1 , p ) ,
        logit(p) <- a[actor] + (bp + bpC*condition)*prosoc_left ,
        a[actor] ~ dnorm(0,10),
        bp ~ dnorm(0,10),
        bpC ~ dnorm(0,10)
    ) ,
    data=d2 , chains=2 , iter=2500 , warmup=500 )

## R code 10.15
unique( d$actor )

## R code 10.16
precis( m10.4 , depth=2 )

## R code 10.17
post <- extract.samples( m10.4 )
str( post )

## R code 10.18
dens( post$a[,2] )

## R code 10.19
chimp <- 1
d.pred <- list(
    pulled_left = rep( 0 , 4 ), # empty outcome
    prosoc_left = c(0,1,0,1),   # right/left/right/left
    condition = c(0,0,1,1),     # control/control/partner/partner
    actor = rep(chimp,4)
)
link.m10.4 <- link( m10.4 , data=d.pred )
pred.p <- apply( link.m10.4 , 2 , mean )
pred.p.PI <- apply( link.m10.4 , 2 , PI )

plot( 0 , 0 , type="n" , xlab="prosoc_left/condition" ,
    ylab="proportion pulled left" , ylim=c(0,1) , xaxt="n" ,
    xlim=c(1,4) , yaxp=c(0,1,2) )
axis( 1 , at=1:4 , labels=c("0/0","1/0","0/1","1/1") )
mtext( paste( "actor" , chimp ) )

p <- by( d$pulled_left ,
    list(d$prosoc_left,d$condition,d$actor) , mean )
lines( 1:4 , as.vector(p[,,chimp]) , col=rangi2 , lwd=2 )

lines( 1:4 , pred.p )
shade( pred.p.PI , 1:4 )

## R code 10.20
data(chimpanzees)
d <- chimpanzees
d.aggregated <- aggregate( d$pulled_left ,
    list(prosoc_left=d$prosoc_left,condition=d$condition,actor=d$actor) ,
    sum )

## R code 10.21
m10.5 <- map(
    alist(
        x ~ dbinom( 18 , p ) ,
        logit(p) <- a + (bp + bpC*condition)*prosoc_left ,
        a ~ dnorm(0,10) ,
        bp ~ dnorm(0,10) ,
        bpC ~ dnorm(0,10)
    ) ,
    data=d.aggregated )

## R code 10.22
library(rethinking)
data(UCBadmit)
d <- UCBadmit

## R code 10.23
d$male <- ifelse( d$applicant.gender=="male" , 1 , 0 )
m10.6 <- map(
    alist(
         admit ~ dbinom( applications , p ) ,
         logit(p) <- a + bm*male ,
         a ~ dnorm(0,10) ,
         bm ~ dnorm(0,10)
    ) ,
    data=d )
m10.7 <- map(
    alist(
         admit ~ dbinom( applications , p ) ,
         logit(p) <- a ,
         a ~ dnorm(0,10)
    ) ,
    data=d )

## R code 10.24
compare( m10.6 , m10.7 )

## R code 10.25
precis(m10.6)

## R code 10.26
post <- extract.samples( m10.6 )
p.admit.male <- logistic( post$a + post$bm )
p.admit.female <- logistic( post$a )
diff.admit <- p.admit.male - p.admit.female
quantile( diff.admit , c(0.025,0.5,0.975) )

## R code 10.27
postcheck( m10.6 , n=1e4 )
# draw lines connecting points from same dept
for ( i in 1:6 ) {
    x <- 1 + 2*(i-1)
    y1 <- d$admit[x]/d$applications[x]
    y2 <- d$admit[x+1]/d$applications[x+1]
    lines( c(x,x+1) , c(y1,y2) , col=rangi2 , lwd=2 )
    text( x+0.5 , (y1+y2)/2 + 0.05 , d$dept[x] , cex=0.8 , col=rangi2 )
}

## R code 10.28
# make index
d$dept_id <- coerce_index( d$dept )

# model with unique intercept for each dept
m10.8 <- map(
    alist(
        admit ~ dbinom( applications , p ) ,
        logit(p) <- a[dept_id] ,
        a[dept_id] ~ dnorm(0,10)
    ) , data=d )

# model with male difference as well
m10.9 <- map(
    alist(
        admit ~ dbinom( applications , p ) ,
        logit(p) <- a[dept_id] + bm*male ,
        a[dept_id] ~ dnorm(0,10) ,
        bm ~ dnorm(0,10)
    ) , data=d )

## R code 10.29
compare( m10.6 , m10.7 , m10.8 , m10.9 )

## R code 10.30
precis( m10.9 , depth=2 )

## R code 10.31
m10.9stan <- map2stan( m10.9 , chains=2 , iter=2500 , warmup=500 )
precis(m10.9stan,depth=2)

## R code 10.32
m10.7glm <- glm( cbind(admit,reject) ~ 1 , data=d , family=binomial )
m10.6glm <- glm( cbind(admit,reject) ~ male , data=d , family=binomial )
m10.8glm <- glm( cbind(admit,reject) ~ dept , data=d , family=binomial )
m10.9glm <- glm( cbind(admit,reject) ~ male + dept , data=d ,
    family=binomial )

## R code 10.33
data(chimpanzees)
m10.4glm <- glm(
    pulled_left ~ as.factor(actor) + prosoc_left * condition - condition ,
    data=chimpanzees , family=binomial )

## R code 10.34
glimmer( pulled_left ~ prosoc_left * condition - condition ,
    data=chimpanzees , family=binomial )

## R code 10.35
# outcome and predictor almost perfectly associated
y <- c( rep(0,10) , rep(1,10) )
x <- c( rep(-1,9) , rep(1,11) )
# fit binomial GLM
m.bad <- glm( y ~ x , data=list(y=y,x=x) , family=binomial )
precis(m.bad)

## R code 10.36
m.good <- map(
    alist(
        y ~ dbinom( 1 , p ),
        logit(p) <- a + b*x,
        c(a,b) ~ dnorm(0,10)
    ) , data=list(y=y,x=x) )
precis(m.good)

## R code 10.37
m.good.stan <- map2stan( m.good )
pairs(m.good.stan)

## R code 10.38
y <- rbinom(1e5,1000,1/1000)
c( mean(y) , var(y) )

## R code 10.39
library(rethinking)
data(Kline)
d <- Kline
d

## R code 10.40
d$log_pop <- log(d$population)
d$contact_high <- ifelse( d$contact=="high" , 1 , 0 )

## R code 10.41
m10.10 <- map(
    alist(
        total_tools ~ dpois( lambda ),
        log(lambda) <- a + bp*log_pop +
            bc*contact_high + bpc*contact_high*log_pop,
        a ~ dnorm(0,100),
        c(bp,bc,bpc) ~ dnorm(0,1)
    ),
    data=d )

## R code 10.42
precis(m10.10,corr=TRUE)
plot(precis(m10.10))

## R code 10.43
post <- extract.samples(m10.10)
lambda_high <- exp( post$a + post$bc + (post$bp + post$bpc)*8 )
lambda_low <- exp( post$a + post$bp*8 )

## R code 10.44
diff <- lambda_high - lambda_low
sum(diff > 0)/length(diff)

## R code 10.45
# no interaction
m10.11 <- map(
    alist(
        total_tools ~ dpois( lambda ),
        log(lambda) <- a + bp*log_pop + bc*contact_high,
        a ~ dnorm(0,100),
        c(bp,bc) ~ dnorm( 0 , 1 )
    ), data=d )

## R code 10.46
# no contact rate
m10.12 <- map(
    alist(
        total_tools ~ dpois( lambda ),
        log(lambda) <- a + bp*log_pop,
        a ~ dnorm(0,100),
        bp ~ dnorm( 0 , 1 )
    ), data=d )

# no log-population
m10.13 <- map(
    alist(
        total_tools ~ dpois( lambda ),
        log(lambda) <- a + bc*contact_high,
        a ~ dnorm(0,100),
        bc ~ dnorm( 0 , 1 )
    ), data=d )

## R code 10.47
# intercept only
m10.14 <- map(
    alist(
        total_tools ~ dpois( lambda ),
        log(lambda) <- a,
        a ~ dnorm(0,100)
    ), data=d )

# compare all using WAIC
# adding n=1e4 for more stable WAIC estimates
# will also plot the comparison
( islands.compare <- compare(m10.10,m10.11,m10.12,m10.13,m10.14,n=1e4) )
plot(islands.compare)

## R code 10.48
# make plot of raw data to begin
# point character (pch) indicates contact rate
pch <- ifelse( d$contact_high==1 , 16 , 1 )
plot( d$log_pop , d$total_tools , col=rangi2 , pch=pch ,
    xlab="log-population" , ylab="total tools" )

# sequence of log-population sizes to compute over
log_pop.seq <- seq( from=6 , to=13 , length.out=30 )

# compute trend for high contact islands
d.pred <- data.frame(
    log_pop = log_pop.seq,
    contact_high = 1
)
lambda.pred.h <- ensemble( m10.10 , m10.11 , m10.12 , data=d.pred )
lambda.med <- apply( lambda.pred.h$link , 2 , median )
lambda.PI <- apply( lambda.pred.h$link , 2 , PI )

# plot predicted trend for high contact islands
lines( log_pop.seq , lambda.med , col=rangi2 )
shade( lambda.PI , log_pop.seq , col=col.alpha(rangi2,0.2) )

# compute trend for low contact islands
d.pred <- data.frame(
    log_pop = log_pop.seq,
    contact_high = 0
)
lambda.pred.l <- ensemble( m10.10 , m10.11 , m10.12 , data=d.pred )
lambda.med <- apply( lambda.pred.l$link , 2 , median )
lambda.PI <- apply( lambda.pred.l$link , 2 , PI )

# plot again
lines( log_pop.seq , lambda.med , lty=2 )
shade( lambda.PI , log_pop.seq , col=col.alpha("black",0.1) )

## R code 10.49
m10.10stan <- map2stan( m10.10 , iter=3000 , warmup=1000 , chains=4 )
precis(m10.10stan)

## R code 10.50
# construct centered predictor
d$log_pop_c <- d$log_pop - mean(d$log_pop)

# re-estimate
m10.10stan.c <- map2stan(
    alist(
        total_tools ~ dpois( lambda ) ,
        log(lambda) <- a + bp*log_pop_c + bc*contact_high +
            bcp*log_pop_c*contact_high ,
        a ~ dnorm(0,10) ,
        bp ~ dnorm(0,1) ,
        bc ~ dnorm(0,1) ,
        bcp ~ dnorm(0,1)
    ) ,
    data=d , iter=3000 , warmup=1000 , chains=4 )
precis(m10.10stan.c)

## R code 10.51
num_days <- 30
y <- rpois( num_days , 1.5 )

## R code 10.52
num_weeks <- 4
y_new <- rpois( num_weeks , 0.5*7 )

## R code 10.53
y_all <- c( y , y_new )
exposure <- c( rep(1,30) , rep(7,4) )
monastery <- c( rep(0,30) , rep(1,4) )
d <- data.frame( y=y_all , days=exposure , monastery=monastery )

## R code 10.54
# compute the offset
d$log_days <- log( d$days )

# fit the model
m10.15 <- map(
    alist(
        y ~ dpois( lambda ),
        log(lambda) <- log_days + a + b*monastery,
        a ~ dnorm(0,100),
        b ~ dnorm(0,1)
    ),
    data=d )

## R code 10.55
post <- extract.samples( m10.15 )
lambda_old <- exp( post$a )
lambda_new <- exp( post$a + post$b )
precis( data.frame( lambda_old , lambda_new ) )

## R code 10.56
# simulate career choices among 500 individuals
N <- 500             # number of individuals
income <- 1:3        # expected income of each career
score <- 0.5*income  # scores for each career, based on income
# next line converts scores to probabilities
p <- softmax(score[1],score[2],score[3])

# now simulate choice
# outcome career holds event type values, not counts
career <- rep(NA,N)  # empty vector of choices for each individual
# sample chosen career for each individual
for ( i in 1:N ) career[i] <- sample( 1:3 , size=1 , prob=p )

## R code 10.57
# fit the model, using dcategorical and softmax link
m10.16 <- map(
    alist(
        career ~ dcategorical( softmax(0,s2,s3) ),
        s2 <- b*2,    # linear model for event type 2
        s3 <- b*3,    # linear model for event type 3
        b ~ dnorm(0,5)
    ) ,
    data=list(career=career) )

## R code 10.58
N <- 100
# simulate family incomes for each individual
family_income <- runif(N)
# assign a unique coefficient for each type of event
b <- (1:-1)
career <- rep(NA,N)  # empty vector of choices for each individual
for ( i in 1:N ) {
    score <- 0.5*(1:3) + b*family_income[i]
    p <- softmax(score[1],score[2],score[3])
    career[i] <- sample( 1:3 , size=1 , prob=p )
}

m10.17 <- map(
    alist(
        career ~ dcategorical( softmax(0,s2,s3) ),
        s2 <- a2 + b2*family_income,
        s3 <- a3 + b3*family_income,
        c(a2,a3,b2,b3) ~ dnorm(0,5)
    ) ,
    data=list(career=career,family_income=family_income) )

## R code 10.59
library(rethinking)
data(UCBadmit)
d <- UCBadmit

## R code 10.60
# binomial model of overall admission probability
m_binom <- map(
    alist(
        admit ~ dbinom(applications,p),
        logit(p) <- a,
        a ~ dnorm(0,100)
    ),
    data=d )

# Poisson model of overall admission rate and rejection rate
d$rej <- d$reject # 'reject' is a reserved word
m_pois <- map2stan(
    alist(
        admit ~ dpois(lambda1),
        rej ~ dpois(lambda2),
        log(lambda1) <- a1,
        log(lambda2) <- a2,
        c(a1,a2) ~ dnorm(0,100)
    ),
    data=d , chains=3 , cores=3 )

## R code 10.61
logistic(coef(m_binom))

## R code 10.62
k <- as.numeric(coef(m_pois))
exp(k[1])/(exp(k[1])+exp(k[2]))

## R code 10.63
# simulate
N <- 100
x <- runif(N)
y <- rgeom( N , prob=logistic( -1 + 2*x ) )

# estimate
m10.18 <- map(
    alist(
        y ~ dgeom( p ),
        logit(p) <- a + b*x,
        a ~ dnorm(0,10),
        b ~ dnorm(0,1)
    ),
    data=list(y=y,x=x) )
precis(m10.18)

## R code 11.1
library(rethinking)
data(Trolley)
d <- Trolley

## R code 11.2
simplehist( d$response , xlim=c(1,7) , xlab="response" )

## R code 11.3
# discrete proportion of each response value
pr_k <- table( d$response ) / nrow(d)

# cumsum converts to cumulative proportions
cum_pr_k <- cumsum( pr_k )

# plot
plot( 1:7 , cum_pr_k , type="b" , xlab="response" ,
ylab="cumulative proportion" , ylim=c(0,1) )

## R code 11.4
logit <- function(x) log(x/(1-x)) # convenience function
( lco <- logit( cum_pr_k ) )

## R code 11.5
m11.1 <- map(
    alist(
        response ~ dordlogit( phi , c(a1,a2,a3,a4,a5,a6) ),
        phi <- 0,
        c(a1,a2,a3,a4,a5,a6) ~ dnorm(0,10)
    ) ,
    data=d ,
    start=list(a1=-2,a2=-1,a3=0,a4=1,a5=2,a6=2.5) )

## R code 11.6
precis(m11.1)

## R code 11.7
logistic(coef(m11.1))

## R code 11.8
# note that data with name 'case' not allowed in Stan
# so will pass pruned data list
m11.1stan <- map2stan(
    alist(
        response ~ dordlogit( phi , cutpoints ),
        phi <- 0,
        cutpoints ~ dnorm(0,10)
    ) ,
    data=list(response=d$response),
    start=list(cutpoints=c(-2,-1,0,1,2,2.5)) ,
    chains=2 , cores=2 )

# need depth=2 to show vector of parameters
precis(m11.1stan,depth=2)

## R code 11.9
( pk <- dordlogit( 1:7 , 0 , coef(m11.1) ) )

## R code 11.10
sum( pk*(1:7) )

## R code 11.11
( pk <- dordlogit( 1:7 , 0 , coef(m11.1)-0.5 ) )

## R code 11.12
sum( pk*(1:7) )

## R code 11.13
m11.2 <- map(
    alist(
        response ~ dordlogit( phi , c(a1,a2,a3,a4,a5,a6) ) ,
        phi <- bA*action + bI*intention + bC*contact,
        c(bA,bI,bC) ~ dnorm(0,10),
        c(a1,a2,a3,a4,a5,a6) ~ dnorm(0,10)
    ) ,
    data=d ,
    start=list(a1=-1.9,a2=-1.2,a3=-0.7,a4=0.2,a5=0.9,a6=1.8) )

## R code 11.14
m11.3 <- map(
    alist(
        response ~ dordlogit( phi , c(a1,a2,a3,a4,a5,a6) ) ,
        phi <- bA*action + bI*intention + bC*contact +
            bAI*action*intention + bCI*contact*intention ,
        c(bA,bI,bC,bAI,bCI) ~ dnorm(0,10),
        c(a1,a2,a3,a4,a5,a6) ~ dnorm(0,10)
    ) ,
    data=d ,
    start=list(a1=-1.9,a2=-1.2,a3=-0.7,a4=0.2,a5=0.9,a6=1.8) )

## R code 11.15
coeftab(m11.1,m11.2,m11.3)

## R code 11.16
compare( m11.1 , m11.2 , m11.3 , refresh=0.1 )

## R code 11.17
post <- extract.samples( m11.3 )

## R code 11.18
plot( 1 , 1 , type="n" , xlab="intention" , ylab="probability" ,
    xlim=c(0,1) , ylim=c(0,1) , xaxp=c(0,1,1) , yaxp=c(0,1,2) )

## R code 11.19
kA <- 0     # value for action
kC <- 1     # value for contact
kI <- 0:1   # values of intention to calculate over
for ( s in 1:100 ) {
    p <- post[s,]
    ak <- as.numeric(p[1:6])
    phi <- p$bA*kA + p$bI*kI + p$bC*kC +
        p$bAI*kA*kI + p$bCI*kC*kI
    pk <- pordlogit( 1:6 , a=ak , phi=phi )
    for ( i in 1:6 )
        lines( kI , pk[,i] , col=col.alpha(rangi2,0.1) )
}
mtext( concat( "action=",kA,", contact=",kC ) )

## R code 11.20
# define parameters
prob_drink <- 0.2 # 20% of days
rate_work <- 1    # average 1 manuscript per day

# sample one year of production
N <- 365

# simulate days monks drink
drink <- rbinom( N , 1 , prob_drink )

# simulate manuscripts completed
y <- (1-drink)*rpois( N , rate_work )

## R code 11.21
simplehist( y , xlab="manuscripts completed" , lwd=4 )
zeros_drink <- sum(drink)
zeros_work <- sum(y==0 & drink==0)
zeros_total <- sum(y==0)
lines( c(0,0) , c(zeros_work,zeros_total) , lwd=4 , col=rangi2 )

## R code 11.22
m11.4 <- map(
    alist(
        y ~ dzipois( p , lambda ),
        logit(p) <- ap,
        log(lambda) <- al,
        ap ~ dnorm(0,1),
        al ~ dnorm(0,10)
    ) ,
    data=list(y=y) )
precis(m11.4)

## R code 11.23
logistic(-1.39) # probability drink
exp(0.05)       # rate finish manuscripts, when not drinking

## R code 11.24
dzip <- function( x , p , lambda , log=TRUE ) {
    ll <- ifelse(
        x==0 ,
        p + (1-p)*exp(-lambda) ,
        (1-p)*dpois(x,lambda,FALSE)
    )
    if ( log==TRUE ) ll <- log(ll)
    return(ll)
}

## R code 11.25
pbar <- 0.5
theta <- 5
curve( dbeta2(x,pbar,theta) , from=0 , to=1 ,
    xlab="probability" , ylab="Density" )

## R code 11.26
library(rethinking)
data(UCBadmit)
d <- UCBadmit
m11.5 <- map2stan(
    alist(
        admit ~ dbetabinom(applications,pbar,theta),
        logit(pbar) <- a,
        a ~ dnorm(0,2),
        theta ~ dexp(1)
    ),
    data=d,
    constraints=list(theta="lower=0"),
    start=list(theta=3),
    iter=4000 , warmup=1000 , chains=2 , cores=2 )

## R code 11.27
precis(m11.5)

## R code 11.28
post <- extract.samples(m11.5)
quantile( logistic(post$a) , c(0.025,0.5,0.975) )

## R code 11.29
post <- extract.samples(m11.5)

# draw posterior mean beta distribution
curve( dbeta2(x,mean(logistic(post$a)),mean(post$theta)) , from=0 , to=1 ,
    ylab="Density" , xlab="probability admit", ylim=c(0,3) , lwd=2 )

# draw 100 beta distributions sampled from posterior
for ( i in 1:100 ) {
    p <- logistic( post$a[i] )
    theta <- post$theta[i]
    curve( dbeta2(x,p,theta) , add=TRUE , col=col.alpha("black",0.2) )
}

## R code 11.30
postcheck(m11.5)

## R code 11.31
mu <- 3
theta <- 1
curve( dgamma2(x,mu,theta) , from=0 , to=10 )

## R code 11.32
library(rethinking)
data(Hurricanes)

## R code 12.1
library(rethinking)
data(reedfrogs)
d <- reedfrogs
str(d)

## R code 12.2
library(rethinking)
data(reedfrogs)
d <- reedfrogs

# make the tank cluster variable
d$tank <- 1:nrow(d)

# fit
m12.1 <- map2stan(
    alist(
        surv ~ dbinom( density , p ) ,
        logit(p) <- a_tank[tank] ,
        a_tank[tank] ~ dnorm( 0 , 5 )
    ),
    data=d )

## R code 12.3
m12.2 <- map2stan(
    alist(
        surv ~ dbinom( density , p ) ,
        logit(p) <- a_tank[tank] ,
        a_tank[tank] ~ dnorm( a , sigma ) ,
        a ~ dnorm(0,1) ,
        sigma ~ dcauchy(0,1)
    ), data=d , iter=4000 , chains=4 )

## R code 12.4
compare( m12.1 , m12.2 )

## R code 12.5
# extract Stan samples
post <- extract.samples(m12.2)

# compute median intercept for each tank
# also transform to probability with logistic
d$propsurv.est <- logistic( apply( post$a_tank , 2 , median ) )

# display raw proportions surviving in each tank
plot( d$propsurv , ylim=c(0,1) , pch=16 , xaxt="n" ,
    xlab="tank" , ylab="proportion survival" , col=rangi2 )
axis( 1 , at=c(1,16,32,48) , labels=c(1,16,32,48) )

# overlay posterior medians
points( d$propsurv.est )

# mark posterior median probability across tanks
abline( h=logistic(median(post$a)) , lty=2 )

# draw vertical dividers between tank densities
abline( v=16.5 , lwd=0.5 )
abline( v=32.5 , lwd=0.5 )
text( 8 , 0 , "small tanks" )
text( 16+8 , 0 , "medium tanks" )
text( 32+8 , 0 , "large tanks" )

## R code 12.6
# show first 100 populations in the posterior
plot( NULL , xlim=c(-3,4) , ylim=c(0,0.35) ,
    xlab="log-odds survive" , ylab="Density" )
for ( i in 1:100 )
    curve( dnorm(x,post$a[i],post$sigma[i]) , add=TRUE ,
    col=col.alpha("black",0.2) )

# sample 8000 imaginary tanks from the posterior distribution
sim_tanks <- rnorm( 8000 , post$a , post$sigma )

# transform to probability and visualize
dens( logistic(sim_tanks) , xlab="probability survive" )

## R code 12.7
a <- 1.4
sigma <- 1.5
nponds <- 60
ni <- as.integer( rep( c(5,10,25,35) , each=15 ) )

## R code 12.8
a_pond <- rnorm( nponds , mean=a , sd=sigma )

## R code 12.9
dsim <- data.frame( pond=1:nponds , ni=ni , true_a=a_pond )

## R code 12.10
class(1:3)
class(c(1,2,3))

## R code 12.11
dsim$si <- rbinom( nponds , prob=logistic(dsim$true_a) , size=dsim$ni )

## R code 12.12
dsim$p_nopool <- dsim$si / dsim$ni

## R code 12.13
m12.3 <- map2stan(
    alist(
        si ~ dbinom( ni , p ),
        logit(p) <- a_pond[pond],
        a_pond[pond] ~ dnorm( a , sigma ),
        a ~ dnorm(0,1),
        sigma ~ dcauchy(0,1)
    ),
    data=dsim , iter=1e4 , warmup=1000 )

## R code 12.14
precis(m12.3,depth=2)

## R code 12.15
estimated.a_pond <- as.numeric( coef(m12.3)[1:60] )
dsim$p_partpool <- logistic( estimated.a_pond )

## R code 12.16
dsim$p_true <- logistic( dsim$true_a )

## R code 12.17
nopool_error <- abs( dsim$p_nopool - dsim$p_true )
partpool_error <- abs( dsim$p_partpool - dsim$p_true )

## R code 12.18
plot( 1:60 , nopool_error , xlab="pond" , ylab="absolute error" ,
    col=rangi2 , pch=16 )
points( 1:60 , partpool_error )

## R code 12.19
a <- 1.4
sigma <- 1.5
nponds <- 60
ni <- as.integer( rep( c(5,10,25,35) , each=15 ) )
a_pond <- rnorm( nponds , mean=a , sd=sigma )
dsim <- data.frame( pond=1:nponds , ni=ni , true_a=a_pond )
dsim$si <- rbinom( nponds,prob=logistic( dsim$true_a ),size=dsim$ni )
dsim$p_nopool <- dsim$si / dsim$ni
newdat <- list(si=dsim$si,ni=dsim$ni,pond=1:nponds)
m12.3new <- map2stan( m12.3 , data=newdat , iter=1e4 , warmup=1000 )

## R code 12.20
y1 <- rnorm( 1e4 , 10 , 1 )
y2 <- 10 + rnorm( 1e4 , 0 , 1 )

## R code 12.21
library(rethinking)
data(chimpanzees)
d <- chimpanzees
d$recipient <- NULL     # get rid of NAs

m12.4 <- map2stan(
    alist(
        pulled_left ~ dbinom( 1 , p ) ,
        logit(p) <- a + a_actor[actor] + (bp + bpC*condition)*prosoc_left ,
        a_actor[actor] ~ dnorm( 0 , sigma_actor ),
        a ~ dnorm(0,10),
        bp ~ dnorm(0,10),
        bpC ~ dnorm(0,10),
        sigma_actor ~ dcauchy(0,1)
    ) ,
    data=d , warmup=1000 , iter=5000 , chains=4 , cores=3 )

## R code 12.22
post <- extract.samples(m12.4)
total_a_actor <- sapply( 1:7 , function(actor) post$a + post$a_actor[,actor] )
round( apply(total_a_actor,2,mean) , 2 )

## R code 12.23
# prep data
d$block_id <- d$block  # name 'block' is reserved by Stan

m12.5 <- map2stan(
    alist(
        pulled_left ~ dbinom( 1 , p ),
        logit(p) <- a + a_actor[actor] + a_block[block_id] +
                    (bp + bpc*condition)*prosoc_left,
        a_actor[actor] ~ dnorm( 0 , sigma_actor ),
        a_block[block_id] ~ dnorm( 0 , sigma_block ),
        c(a,bp,bpc) ~ dnorm(0,10),
        sigma_actor ~ dcauchy(0,1),
        sigma_block ~ dcauchy(0,1)
    ) ,
    data=d, warmup=1000 , iter=6000 , chains=4 , cores=3 )

## R code 12.24
precis(m12.5,depth=2) # depth=2 displays varying effects
plot(precis(m12.5,depth=2)) # also plot

## R code 12.25
post <- extract.samples(m12.5)
dens( post$sigma_block , xlab="sigma" , xlim=c(0,4) )
dens( post$sigma_actor , col=rangi2 , lwd=2 , add=TRUE )
text( 2 , 0.85 , "actor" , col=rangi2 )
text( 0.75 , 2 , "block" )

## R code 12.26
compare(m12.4,m12.5)

## R code 12.27
chimp <- 2
d.pred <- list(
    prosoc_left = c(0,1,0,1),   # right/left/right/left
    condition = c(0,0,1,1),     # control/control/partner/partner
    actor = rep(chimp,4)
)
link.m12.4 <- link( m12.4 , data=d.pred )
pred.p <- apply( link.m12.4 , 2 , mean )
pred.p.PI <- apply( link.m12.4 , 2 , PI )

## R code 12.28
post <- extract.samples(m12.4)
str(post)

## R code 12.29
dens( post$a_actor[,5] )

## R code 12.30
p.link <- function( prosoc_left , condition , actor ) {
    logodds <- with( post ,
        a + a_actor[,actor] + (bp + bpC * condition) * prosoc_left
    )
    return( logistic(logodds) )
}

## R code 12.31
prosoc_left <- c(0,1,0,1)
condition <- c(0,0,1,1)
pred.raw <- sapply( 1:4 , function(i) p.link(prosoc_left[i],condition[i],2) )
pred.p <- apply( pred.raw , 2 , mean )
pred.p.PI <- apply( pred.raw , 2 , PI )

## R code 12.32
d.pred <- list(
    prosoc_left = c(0,1,0,1),   # right/left/right/left
    condition = c(0,0,1,1),     # control/control/partner/partner
    actor = rep(2,4) )          # placeholder

## R code 12.33
# replace varying intercept samples with zeros
# 1000 samples by 7 actors
a_actor_zeros <- matrix(0,1000,7)

## R code 12.34
# fire up link
# note use of replace list
link.m12.4 <- link( m12.4 , n=1000 , data=d.pred ,
    replace=list(a_actor=a_actor_zeros) )

# summarize and plot
pred.p.mean <- apply( link.m12.4 , 2 , mean )
pred.p.PI <- apply( link.m12.4 , 2 , PI , prob=0.8 )
plot( 0 , 0 , type="n" , xlab="prosoc_left/condition" ,
    ylab="proportion pulled left" , ylim=c(0,1) , xaxt="n" ,
    xlim=c(1,4) )
axis( 1 , at=1:4 , labels=c("0/0","1/0","0/1","1/1") )
lines( 1:4 , pred.p.mean )
shade( pred.p.PI , 1:4 )

## R code 12.35
# replace varying intercept samples with simulations
post <- extract.samples(m12.4)
a_actor_sims <- rnorm(7000,0,post$sigma_actor)
a_actor_sims <- matrix(a_actor_sims,1000,7)

## R code 12.36
link.m12.4 <- link( m12.4 , n=1000 , data=d.pred ,
    replace=list(a_actor=a_actor_sims) )

## R code 12.37
post <- extract.samples(m12.4)
sim.actor <- function(i) {
    sim_a_actor <- rnorm( 1 , 0 , post$sigma_actor[i] )
    P <- c(0,1,0,1)
    C <- c(0,0,1,1)
    p <- logistic(
        post$a[i] +
        sim_a_actor +
        (post$bp[i] + post$bpC[i]*C)*P
    )
    return(p)
}

## R code 12.38
# empty plot
plot( 0 , 0 , type="n" , xlab="prosoc_left/condition" ,
    ylab="proportion pulled left" , ylim=c(0,1) , xaxt="n" , xlim=c(1,4) )
axis( 1 , at=1:4 , labels=c("0/0","1/0","0/1","1/1") )

# plot 50 simulated actors
for ( i in 1:50 ) lines( 1:4 , sim.actor(i) , col=col.alpha("black",0.5) )

## R code 12.39
# prep data
library(rethinking)
data(Kline)
d <- Kline
d$logpop <- log(d$population)
d$society <- 1:10

# fit model
m12.6 <- map2stan(
    alist(
        total_tools ~ dpois(mu),
        log(mu) <- a + a_society[society] + bp*logpop,
        a ~ dnorm(0,10),
        bp ~ dnorm(0,1),
        a_society[society] ~ dnorm(0,sigma_society),
        sigma_society ~ dcauchy(0,1)
    ),
    data=d ,
    iter=4000 , chains=3 )

## R code 12.40
post <- extract.samples(m12.6)
d.pred <- list(
    logpop = seq(from=6,to=14,length.out=30),
    society = rep(1,30)
)
a_society_sims <- rnorm(20000,0,post$sigma_society)
a_society_sims <- matrix(a_society_sims,2000,10)
link.m12.6 <- link( m12.6 , n=2000 , data=d.pred ,
    replace=list(a_society=a_society_sims) )

## R code 12.41
# plot raw data
plot( d$logpop , d$total_tools , col=rangi2 , pch=16 ,
    xlab="log population" , ylab="total tools" )

# plot posterior median
mu.median <- apply( link.m12.6 , 2 , median )
lines( d.pred$logpop , mu.median )

# plot 97%, 89%, and 67% intervals (all prime numbers)
mu.PI <- apply( link.m12.6 , 2 , PI , prob=0.97 )
shade( mu.PI , d.pred$logpop )
mu.PI <- apply( link.m12.6 , 2 , PI , prob=0.89 )
shade( mu.PI , d.pred$logpop )
mu.PI <- apply( link.m12.6 , 2 , PI , prob=0.67 )
shade( mu.PI , d.pred$logpop )

## R code 12.42
sort(unique(d$district))

## R code 12.43
d$district_id <- as.integer(as.factor(d$district))
sort(unique(d$district_id))

## R code 13.1
a <- 3.5            # average morning wait time
b <- (-1)           # average difference afternoon wait time
sigma_a <- 1        # std dev in intercepts
sigma_b <- 0.5      # std dev in slopes
rho <- (-0.7)       # correlation between intercepts and slopes

## R code 13.2
Mu <- c( a , b )

## R code 13.3
cov_ab <- sigma_a*sigma_b*rho
Sigma <- matrix( c(sigma_a^2,cov_ab,cov_ab,sigma_b^2) , ncol=2 )

## R code 13.4
matrix( c(1,2,3,4) , nrow=2 , ncol=2 )

## R code 13.5
sigmas <- c(sigma_a,sigma_b) # standard deviations
Rho <- matrix( c(1,rho,rho,1) , nrow=2 ) # correlation matrix

# now matrix multiply to get covariance matrix
Sigma <- diag(sigmas) %*% Rho %*% diag(sigmas)

## R code 13.6
N_cafes <- 20

## R code 13.7
library(MASS)
set.seed(5) # used to replicate example
vary_effects <- mvrnorm( N_cafes , Mu , Sigma )

## R code 13.8
a_cafe <- vary_effects[,1]
b_cafe <- vary_effects[,2]

## R code 13.9
plot( a_cafe , b_cafe , col=rangi2 ,
    xlab="intercepts (a_cafe)" , ylab="slopes (b_cafe)" )

# overlay population distribution
library(ellipse)
for ( l in c(0.1,0.3,0.5,0.8,0.99) )
    lines(ellipse(Sigma,centre=Mu,level=l),col=col.alpha("black",0.2))

## R code 13.10
N_visits <- 10
afternoon <- rep(0:1,N_visits*N_cafes/2)
cafe_id <- rep( 1:N_cafes , each=N_visits )
mu <- a_cafe[cafe_id] + b_cafe[cafe_id]*afternoon
sigma <- 0.5  # std dev within cafes
wait <- rnorm( N_visits*N_cafes , mu , sigma )
d <- data.frame( cafe=cafe_id , afternoon=afternoon , wait=wait )

## R code 13.11
R <- rlkjcorr( 1e4 , K=2 , eta=2 )
dens( R[,1,2] , xlab="correlation" )

## R code 13.12
m13.1 <- map2stan(
    alist(
        wait ~ dnorm( mu , sigma ),
        mu <- a_cafe[cafe] + b_cafe[cafe]*afternoon,
        c(a_cafe,b_cafe)[cafe] ~ dmvnorm2(c(a,b),sigma_cafe,Rho),
        a ~ dnorm(0,10),
        b ~ dnorm(0,10),
        sigma_cafe ~ dcauchy(0,2),
        sigma ~ dcauchy(0,2),
        Rho ~ dlkjcorr(2)
    ) ,
    data=d ,
    iter=5000 , warmup=2000 , chains=2 )

## R code 13.13
post <- extract.samples(m13.1)
dens( post$Rho[,1,2] )

## R code 13.14
# compute unpooled estimates directly from data
a1 <- sapply( 1:N_cafes ,
        function(i) mean(wait[cafe_id==i & afternoon==0]) )
b1 <- sapply( 1:N_cafes ,
        function(i) mean(wait[cafe_id==i & afternoon==1]) ) - a1

# extract posterior means of partially pooled estimates
post <- extract.samples(m13.1)
a2 <- apply( post$a_cafe , 2 , mean )
b2 <- apply( post$b_cafe , 2 , mean )

# plot both and connect with lines
plot( a1 , b1 , xlab="intercept" , ylab="slope" ,
    pch=16 , col=rangi2 , ylim=c( min(b1)-0.1 , max(b1)+0.1 ) ,
    xlim=c( min(a1)-0.1 , max(a1)+0.1 ) )
points( a2 , b2 , pch=1 )
for ( i in 1:N_cafes ) lines( c(a1[i],a2[i]) , c(b1[i],b2[i]) )

## R code 13.15
# compute posterior mean bivariate Gaussian
Mu_est <- c( mean(post$a) , mean(post$b) )
rho_est <- mean( post$Rho[,1,2] )
sa_est <- mean( post$sigma_cafe[,1] )
sb_est <- mean( post$sigma_cafe[,2] )
cov_ab <- sa_est*sb_est*rho_est
Sigma_est <- matrix( c(sa_est^2,cov_ab,cov_ab,sb_est^2) , ncol=2 )

# draw contours
library(ellipse)
for ( l in c(0.1,0.3,0.5,0.8,0.99) )
    lines(ellipse(Sigma_est,centre=Mu_est,level=l),
        col=col.alpha("black",0.2))

## R code 13.16
# convert varying effects to waiting times
wait_morning_1 <- (a1)
wait_afternoon_1 <- (a1 + b1)
wait_morning_2 <- (a2)
wait_afternoon_2 <- (a2 + b2)

## R code 13.17
library(rethinking)
data(UCBadmit)
d <- UCBadmit
d$male <- ifelse( d$applicant.gender=="male" , 1 , 0 )
d$dept_id <- coerce_index( d$dept )

## R code 13.18
m13.2 <- map2stan(
    alist(
        admit ~ dbinom( applications , p ),
        logit(p) <- a_dept[dept_id] + bm*male,
        a_dept[dept_id] ~ dnorm( a , sigma_dept ),
        a ~ dnorm(0,10),
        bm ~ dnorm(0,1),
        sigma_dept ~ dcauchy(0,2)
    ) ,
    data=d , warmup=500 , iter=4500 , chains=3 )
precis( m13.2 , depth=2 ) # depth=2 to display vector parameters

## R code 13.19
m13.3 <- map2stan(
    alist(
        admit ~ dbinom( applications , p ),
        logit(p) <- a_dept[dept_id] +
                    bm_dept[dept_id]*male,
        c(a_dept,bm_dept)[dept_id] ~ dmvnorm2( c(a,bm) , sigma_dept , Rho ),
        a ~ dnorm(0,10),
        bm ~ dnorm(0,1),
        sigma_dept ~ dcauchy(0,2),
        Rho ~ dlkjcorr(2)
    ) ,
    data=d , warmup=1000 , iter=5000 , chains=4 , cores=3 )

## R code 13.20
plot( precis(m13.3,pars=c("a_dept","bm_dept"),depth=2) )

## R code 13.21
m13.4 <- map2stan(
    alist(
        admit ~ dbinom( applications , p ),
        logit(p) <- a_dept[dept_id],
        a_dept[dept_id] ~ dnorm( a , sigma_dept ),
        a ~ dnorm(0,10),
        sigma_dept ~ dcauchy(0,2)
    ) ,
    data=d , warmup=500 , iter=4500 , chains=3 )

compare( m13.2 , m13.3 , m13.4 )

## R code 13.22
library(rethinking)
data(chimpanzees)
d <- chimpanzees
d$recipient <- NULL
d$block_id <- d$block

m13.6 <- map2stan(
    alist(
        # likeliood
        pulled_left ~ dbinom(1,p),

        # linear models
        logit(p) <- A + (BP + BPC*condition)*prosoc_left,
        A <- a + a_actor[actor] + a_block[block_id],
        BP <- bp + bp_actor[actor] + bp_block[block_id],
        BPC <- bpc + bpc_actor[actor] + bpc_block[block_id],

        # adaptive priors
        c(a_actor,bp_actor,bpc_actor)[actor] ~
                                dmvnorm2(0,sigma_actor,Rho_actor),
        c(a_block,bp_block,bpc_block)[block_id] ~
                                dmvnorm2(0,sigma_block,Rho_block),

        # fixed priors
        c(a,bp,bpc) ~ dnorm(0,1),
        sigma_actor ~ dcauchy(0,2),
        sigma_block ~ dcauchy(0,2),
        Rho_actor ~ dlkjcorr(4),
        Rho_block ~ dlkjcorr(4)
    ) , data=d , iter=5000 , warmup=1000 , chains=3 , cores=3 )

## R code 13.23
m13.6NC <- map2stan(
    alist(
        pulled_left ~ dbinom(1,p),
        logit(p) <- A + (BP + BPC*condition)*prosoc_left,
        A <- a + a_actor[actor] + a_block[block_id],
        BP <- bp + bp_actor[actor] + bp_block[block_id],
        BPC <- bpc + bpc_actor[actor] + bpc_block[block_id],
        # adaptive NON-CENTERED priors
        c(a_actor,bp_actor,bpc_actor)[actor] ~
                                dmvnormNC(sigma_actor,Rho_actor),
        c(a_block,bp_block,bpc_block)[block_id] ~
                                dmvnormNC(sigma_block,Rho_block),
        c(a,bp,bpc) ~ dnorm(0,1),
        sigma_actor ~ dcauchy(0,2),
        sigma_block ~ dcauchy(0,2),
        Rho_actor ~ dlkjcorr(4),
        Rho_block ~ dlkjcorr(4)
    ) , data=d , iter=5000 , warmup=1000 , chains=3 , cores=3 )

## R code 13.24
# extract n_eff values for each model
neff_c <- precis(m13.6,2)@output$n_eff
neff_nc <- precis(m13.6NC,2)@output$n_eff
# plot distributions
boxplot( list( 'm13.6'=neff_c , 'm13.6NC'=neff_nc ) ,
    ylab="effective samples" , xlab="model" )

## R code 13.25
precis( m13.6NC , depth=2 , pars=c("sigma_actor","sigma_block") )

## R code 13.26
p <- link(m13.6NC)
str(p)

## R code 13.27
compare( m13.6NC , m12.5 )

## R code 13.28
m13.6nc1 <- map2stan(
    alist(
        pulled_left ~ dbinom(1,p),

        # linear models
        logit(p) <- A + (BP + BPC*condition)*prosoc_left,
        A <- a + za_actor[actor]*sigma_actor[1] +
                 za_block[block_id]*sigma_block[1],
        BP <- bp + zbp_actor[actor]*sigma_actor[2] +
                   zbp_block[block_id]*sigma_block[2],
        BPC <- bpc + zbpc_actor[actor]*sigma_actor[3] +
                     zbpc_block[block_id]*sigma_block[3],

        # adaptive priors
        c(za_actor,zbp_actor,zbpc_actor)[actor] ~ dmvnorm(0,Rho_actor),
        c(za_block,zbp_block,zbpc_block)[block_id] ~ dmvnorm(0,Rho_block),

        # fixed priors
        c(a,bp,bpc) ~ dnorm(0,1),
        sigma_actor ~ dcauchy(0,2),
        sigma_block ~ dcauchy(0,2),
        Rho_actor ~ dlkjcorr(4),
        Rho_block ~ dlkjcorr(4)
    ) ,
    data=d ,
    start=list( sigma_actor=c(1,1,1), sigma_block=c(1,1,1) ),
    constraints=list( sigma_actor="lower=0", sigma_block="lower=0" ),
    types=list( Rho_actor="corr_matrix", Rho_block="corr_matrix" ),
    iter=5000 , warmup=1000 , chains=3 , cores=3 )

## R code 13.29
# load the distance matrix
library(rethinking)
data(islandsDistMatrix)

# display short column names, so fits on screen
Dmat <- islandsDistMatrix
colnames(Dmat) <- c("Ml","Ti","SC","Ya","Fi","Tr","Ch","Mn","To","Ha")
round(Dmat,1)

## R code 13.30
# linear
curve( exp(-1*x) , from=0 , to=4 , lty=2 ,
    xlab="distance" , ylab="correlation" )

# squared
curve( exp(-1*x^2) , add=TRUE )

## R code 13.31
data(Kline2) # load the ordinary data, now with coordinates
d <- Kline2
d$society <- 1:10 # index observations

m13.7 <- map2stan(
    alist(
        total_tools ~ dpois(lambda),
        log(lambda) <- a + g[society] + bp*logpop,
        g[society] ~ GPL2( Dmat , etasq , rhosq , 0.01 ),
        a ~ dnorm(0,10),
        bp ~ dnorm(0,1),
        etasq ~ dcauchy(0,1),
        rhosq ~ dcauchy(0,1)
    ),
    data=list(
        total_tools=d$total_tools,
        logpop=d$logpop,
        society=d$society,
        Dmat=islandsDistMatrix),
    warmup=2000 , iter=1e4 , chains=4 )

## R code 13.32
precis(m13.7,depth=2)

## R code 13.33
post <- extract.samples(m13.7)

# plot the posterior median covariance function
curve( median(post$etasq)*exp(-median(post$rhosq)*x^2) , from=0 , to=10 ,
    xlab="distance (thousand km)" , ylab="covariance" , ylim=c(0,1) ,
    yaxp=c(0,1,4) , lwd=2 )

# plot 100 functions sampled from posterior
for ( i in 1:100 )
    curve( post$etasq[i]*exp(-post$rhosq[i]*x^2) , add=TRUE ,
        col=col.alpha("black",0.2) )

## R code 13.34
# compute posterior median covariance among societies
K <- matrix(0,nrow=10,ncol=10)
for ( i in 1:10 )
    for ( j in 1:10 )
        K[i,j] <- median(post$etasq) *
                  exp( -median(post$rhosq) * islandsDistMatrix[i,j]^2 )
diag(K) <- median(post$etasq) + 0.01

## R code 13.35
# convert to correlation matrix
Rho <- round( cov2cor(K) , 2 )
# add row/col names for convenience
colnames(Rho) <- c("Ml","Ti","SC","Ya","Fi","Tr","Ch","Mn","To","Ha")
rownames(Rho) <- colnames(Rho)
Rho

## R code 13.36
# scale point size to logpop
psize <- d$logpop / max(d$logpop)
psize <- exp(psize*1.5)-2

# plot raw data and labels
plot( d$lon2 , d$lat , xlab="longitude" , ylab="latitude" ,
    col=rangi2 , cex=psize , pch=16 , xlim=c(-50,30) )
labels <- as.character(d$culture)
text( d$lon2 , d$lat , labels=labels , cex=0.7 , pos=c(2,4,3,3,4,1,3,2,4,2) )

# overlay lines shaded by Rho
for( i in 1:10 )
    for ( j in 1:10 )
        if ( i < j )
            lines( c( d$lon2[i],d$lon2[j] ) , c( d$lat[i],d$lat[j] ) ,
                lwd=2 , col=col.alpha("black",Rho[i,j]^2) )

## R code 13.37
# compute posterior median relationship, ignoring distance
logpop.seq <- seq( from=6 , to=14 , length.out=30 )
lambda <- sapply( logpop.seq , function(lp) exp( post$a + post$bp*lp ) )
lambda.median <- apply( lambda , 2 , median )
lambda.PI80 <- apply( lambda , 2 , PI , prob=0.8 )

# plot raw data and labels
plot( d$logpop , d$total_tools , col=rangi2 , cex=psize , pch=16 ,
    xlab="log population" , ylab="total tools" )
text( d$logpop , d$total_tools , labels=labels , cex=0.7 ,
    pos=c(4,3,4,2,2,1,4,4,4,2) )

# display posterior predictions
lines( logpop.seq , lambda.median , lty=2 )
lines( logpop.seq , lambda.PI80[1,] , lty=2 )
lines( logpop.seq , lambda.PI80[2,] , lty=2 )

# overlay correlations
for( i in 1:10 )
    for ( j in 1:10 )
        if ( i < j )
            lines( c( d$logpop[i],d$logpop[j] ) ,
                   c( d$total_tools[i],d$total_tools[j] ) ,
                   lwd=2 , col=col.alpha("black",Rho[i,j]^2) )

## R code 13.38
S <- matrix( c( sa^2 , sa*sb*rho , sa*sb*rho , sb^2 ) , nrow=2 )

## R code 14.1
# simulate a pancake and return randomly ordered sides
sim_pancake <- function() {
    pancake <- sample(1:3,1)
    sides <- matrix(c(1,1,1,0,0,0),2,3)[,pancake]
    sample(sides)
}

# sim 10,000 pancakes
pancakes <- replicate( 1e4 , sim_pancake() )
up <- pancakes[1,]
down <- pancakes[2,]

# compute proportion 1/1 (BB) out of all 1/1 and 1/0
num_11_10 <- sum( up==1 )
num_11 <- sum( up==1 & down==1 )
num_11/num_11_10

## R code 14.2
library(rethinking)
data(WaffleDivorce)
d <- WaffleDivorce

# points
plot( d$Divorce ~ d$MedianAgeMarriage , ylim=c(4,15) ,
    xlab="Median age marriage" , ylab="Divorce rate" )

# standard errors
for ( i in 1:nrow(d) ) {
    ci <- d$Divorce[i] + c(-1,1)*d$Divorce.SE[i]
    x <- d$MedianAgeMarriage[i]
    lines( c(x,x) , ci )
}

## R code 14.3
dlist <- list(
    div_obs=d$Divorce,
    div_sd=d$Divorce.SE,
    R=d$Marriage,
    A=d$MedianAgeMarriage
)

m14.1 <- map2stan(
    alist(
        div_est ~ dnorm(mu,sigma),
        mu <- a + bA*A + bR*R,
        div_obs ~ dnorm(div_est,div_sd),
        a ~ dnorm(0,10),
        bA ~ dnorm(0,10),
        bR ~ dnorm(0,10),
        sigma ~ dcauchy(0,2.5)
    ) ,
    data=dlist ,
    start=list(div_est=dlist$div_obs) ,
    WAIC=FALSE , iter=5000 , warmup=1000 , chains=2 , cores=2 ,
    control=list(adapt_delta=0.95) )

## R code 14.4
precis( m14.1 , depth=2 )

## R code 14.5
dlist <- list(
    div_obs=d$Divorce,
    div_sd=d$Divorce.SE,
    mar_obs=d$Marriage,
    mar_sd=d$Marriage.SE,
    A=d$MedianAgeMarriage )

m14.2 <- map2stan(
    alist(
        div_est ~ dnorm(mu,sigma),
        mu <- a + bA*A + bR*mar_est[i],
        div_obs ~ dnorm(div_est,div_sd),
        mar_obs ~ dnorm(mar_est,mar_sd),
        a ~ dnorm(0,10),
        bA ~ dnorm(0,10),
        bR ~ dnorm(0,10),
        sigma ~ dcauchy(0,2.5)
    ) ,
    data=dlist ,
    start=list(div_est=dlist$div_obs,mar_est=dlist$mar_obs) ,
    WAIC=FALSE , iter=5000 , warmup=1000 , chains=3 , cores=3 ,
    control=list(adapt_delta=0.95) )

## R code 14.6
library(rethinking)
data(milk)
d <- milk
d$neocortex.prop <- d$neocortex.perc / 100
d$logmass <- log(d$mass)

## R code 14.7
# prep data
data_list <- list(
    kcal = d$kcal.per.g,
    neocortex = d$neocortex.prop,
    logmass = d$logmass )

# fit model
m14.3 <- map2stan(
    alist(
        kcal ~ dnorm(mu,sigma),
        mu <- a + bN*neocortex + bM*logmass,
        neocortex ~ dnorm(nu,sigma_N),
        a ~ dnorm(0,100),
        c(bN,bM) ~ dnorm(0,10),
        nu ~ dnorm(0.5,1),
        sigma_N ~ dcauchy(0,1),
        sigma ~ dcauchy(0,1)
    ) ,
    data=data_list , iter=1e4 , chains=2 )

## R code 14.8
precis(m14.3,depth=2)

## R code 14.9
# prep data
dcc <- d[ complete.cases(d$neocortex.prop) , ]
data_list_cc <- list(
    kcal = dcc$kcal.per.g,
    neocortex = dcc$neocortex.prop,
    logmass = dcc$logmass )

# fit model
m14.3cc <- map2stan(
    alist(
        kcal ~ dnorm(mu,sigma),
        mu <- a + bN*neocortex + bM*logmass,
        a ~ dnorm(0,100),
        c(bN,bM) ~ dnorm(0,10),
        sigma ~ dcauchy(0,1)
    ) ,
    data=data_list_cc , iter=1e4 , chains=2 )
precis(m14.3cc)

## R code 14.10
m14.4 <- map2stan(
    alist(
        kcal ~ dnorm(mu,sigma),
        mu <- a + bN*neocortex + bM*logmass,
        neocortex ~ dnorm(nu,sigma_N),
        nu <- a_N + gM*logmass,
        a ~ dnorm(0,100),
        c(bN,bM,gM) ~ dnorm(0,10),
        a_N ~ dnorm(0.5,1),
        sigma_N ~ dcauchy(0,1),
        sigma ~ dcauchy(0,1)
    ) ,
    data=data_list , iter=1e4 , chains=2 )
precis(m14.4,depth=2)

## R code 14.11
nc_missing <- ifelse( is.na(d$neocortex.prop) , 1 , 0 )
nc_missing <- sapply( 1:length(nc_missing) ,
    function(n) nc_missing[n]*sum(nc_missing[1:n]) )
nc_missing

## R code 14.12
nc <- ifelse( is.na(d$neocortex.prop) , -1 , d$neocortex.prop )

## R code 14.13
model_code <- '
data{
    int N;
    int nc_num_missing;
    vector[N] kcal;
    real neocortex[N];
    vector[N] logmass;
    int nc_missing[N];
}
parameters{
    real alpha;
    real<lower=0> sigma;
    real bN;
    real bM;
    vector[nc_num_missing] nc_impute;
    real mu_nc;
    real<lower=0> sigma_nc;
}
model{
    vector[N] mu;
    vector[N] nc_merged;
    alpha ~ normal(0,10);
    bN ~ normal(0,10);
    bM ~ normal(0,10);
    mu_nc ~ normal(0.5,1);
    sigma ~ cauchy(0,1);
    sigma_nc ~ cauchy(0,1);
    // merge missing and observed
    for ( i in 1:N ) {
        nc_merged[i] <- neocortex[i];
        if ( nc_missing[i] > 0 ) nc_merged[i] <- nc_impute[nc_missing[i]];
    }
    // imputation
    nc_merged ~ normal( mu_nc , sigma_nc );
    // regression
    mu <- alpha + bN*nc_merged + bM*logmass;
    kcal ~ normal( mu , sigma );
}'

## R code 14.14
data_list <- list(
    N = nrow(d),
    kcal = d$kcal.per.g,
    neocortex = nc,
    logmass = d$logmass,
    nc_missing = nc_missing,
    nc_num_missing = max(nc_missing)
)
start <- list(
    alpha=mean(d$kcal.per.g), sigma=sd(d$kcal.per.g),
    bN=0, bM=0, mu_nc=0.68, sigma_nc=0.06,
    nc_impute=rep( 0.5 , max(nc_missing) )
)
library(rstan)
m14.3stan <- stan( model_code=model_code , data=data_list , init=list(start) ,
    iter=1e4 , chains=1 )

## R code 14.15
set.seed(100)
x <- c( rnorm(10) , NA )
y <- c( rnorm(10,x) , 100 )
d <- list(x=x,y=y)

