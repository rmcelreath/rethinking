# distributions

# convert to sigma from tau
tau2sig <- function( tau ) 1/abs(tau)

# Gaussian density that takes tau instead of sigma
dnormtau <- function( x , mean=0 , tau=1 , log=FALSE ) dnorm( x , mean=mean , sd=tau2sig( tau ) , log=log )

rnormtau <- function( n , mean=0 , tau=1 ) rnorm( n=n , mean=mean , sd=tau2sig(tau) )

# ordered categorical density functions
logistic <- function( x ) {
    p <- 1 / ( 1 + exp( -x ) )
    p <- ifelse( x==Inf , 1 , p )
    p
}

pordlogit <- function( x , a , phi , log=FALSE ) {
    a <- c( as.numeric(a) , Inf )
    if ( length(phi) == 1 ) {
        p <- logistic( a[x] - phi )
    } else {
        p <- matrix( NA , ncol=length(x) , nrow=length(phi) )
        for ( i in 1:length(phi) ) {
            p[i,] <- logistic( a[x] - phi[i] )
        }
    }
    if ( log==TRUE ) p <- log(p)
    p
}

dordlogit <- function( x , a , phi , log=FALSE ) {
    a <- c( as.numeric(a) , Inf )
    p <- logistic( a[x] - phi )
    na <- c( -Inf , a )
    np <- logistic( na[x] - phi )
    p <- p - np
    if ( log==TRUE ) p <- log(p)
    p
}

rordlogit <- function( n , a , phi=0 ) {
    a <- c( as.numeric(a) , Inf )
    k <- 1:length(a)
    p <- dordlogit( k , a=a , phi=phi , log=FALSE )
    y <- sample( k , size=n , replace=TRUE , prob=p )
    y
}

# beta density functions parameterized as prob,theta

dbeta2 <- function( x , prob , theta , log=FALSE ) {
    a <- prob * theta
    b <- (1-prob) * theta
    dbeta( x , shape1=a , shape2=b , log=log )
}

rbeta2 <- function( n , prob , theta ) {
    a <- prob * theta
    b <- (1-prob) * theta
    rbeta( n , shape1=a , shape2=b )
}

# based on Ben Bolker's dbetabinom in emdbook package
dbetabinom <- function (x, prob, size, theta, shape1, shape2, log = FALSE) 
{
    if (missing(prob) && !missing(shape1) && !missing(shape2)) {
        prob <- shape1/(shape1 + shape2)
        theta <- shape1 + shape2
    }
    v <- lfactorial(size) - lfactorial(x) - lfactorial(size - 
        x) - lbeta(theta * (1 - prob), theta * prob) + lbeta(size - 
        x + theta * (1 - prob), x + theta * prob)
    if (any(n <- nonint(x))) {
        warning("non-integer x detected; returning zero probability")
        v[n] <- -Inf
    }
    if (log) 
        v
    else exp(v)
}

# gamma-poisson functions

dgamma2 <- function( x , mu , scale , log=FALSE ) {
    dgamma( x , shape=mu / scale , scale=scale , log=log )
}

rgamma2 <- function( n , mu , scale ) {
    rgamma( n , shape=mu/scale , scale=scale )
}

dgampois <- function( x , mu=NULL , shape=NULL , scale=NULL , rate=NULL , log=FALSE ) {
    if ( !is.null(rate) ) scale <- 1/rate
    if ( !is.null(mu) ) shape <- mu / scale
    prob <- 1 / ( 1 + scale )
    dnbinom( x , size=shape , prob=prob , log=log )
}

rgampois <- function( n , mu=NULL , shape=NULL , scale=NULL , rate=NULL ) {
    if ( !is.null(rate) ) scale <- 1/rate
    if ( !is.null(mu) ) shape <- mu / scale
    prob <- 1 / ( 1 + scale )
    rnbinom( n , size=shape , prob=prob )
}

###
# zero-inflated gamma functions
dzerogamma <- function( x , prob , mu , scale , log=TRUE ) {
    p1 <- dgamma2( x , mu=mu , s=scale , log=FALSE )
    p <- ifelse( x==0 , prob , (1-prob)*p1 )
    if ( log==TRUE ) p <- log(p)
    p
}

# laplace (double exponential)
dlaplace <- function(x,lambda,location=0,log=FALSE) {
    # f(y) = (1/(2b)) exp( -|y-a|/b )
    l <- (1/(2*lambda))*exp( -abs(x-location)/lambda )
    if ( log==TRUE ) l <- log(l)
    l
}

rlaplace <- function(n,lambda,location=0) {
    # generate an exponential sample, then generate a random sign, then add location
    y <- rexp(n=n,rate=lambda) * sample( c(-1,1) , size=n , replace=TRUE )
    y <- y + location
    y
}

### softmax rule for multinomial
multilogistic <- function( x , lambda=1 , diff=TRUE , log=FALSE ) {
    if ( diff==TRUE ) x <- x - min(x)
    f <- exp( lambda * x )
    if ( log==FALSE ) result <- f / sum(f)
    if ( log==TRUE ) result <- log(f) - log(sum(f))
    result
}