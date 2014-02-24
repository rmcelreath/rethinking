# distributions

# ordered categorical density functions
logistic <- function( x ) {
    p <- 1 / ( 1 + exp( -x ) )
    p <- ifelse( x==Inf , 1 , p )
    p
}

pordlogit <- function( x , phi , a , log=FALSE ) {
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

dordlogit <- function( x , phi , a , log=FALSE ) {
    a <- c( as.numeric(a) , Inf )
    p <- logistic( a[x] - phi )
    na <- c( -Inf , a )
    np <- logistic( na[x] - phi )
    p <- p - np
    if ( log==TRUE ) p <- log(p)
    p
}

rordlogit <- function( n , phi=0 , a ) {
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

# Ben Bolker's dbetabinom from emdbook package
dbetabinom <- function (x, size, prob, theta, shape1, shape2, log = FALSE) 
{
    if (missing(prob) && !missing(shape1) && !missing(shape2)) {
        prob <- shape1/(shape1 + shape2)
        theta <- shape1 + shape2
    }
    v <- lfactorial(size) - lfactorial(x) - lfactorial(size - 
        x) - lbeta(theta * (1 - prob), theta * prob) + lbeta(size - 
        x + theta * (1 - prob), x + theta * prob)
    nonint <- function(x) (abs((x) - floor((x)+0.5)) > 1e-7)
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

dgampois <- function( x , mu , scale , log=FALSE ) {
    shape <- mu / scale
    prob <- 1 / ( 1 + scale )
    dnbinom( x , size=shape , prob=prob , log=log )
}

rgampois <- function( n , mu , scale ) {
    shape <- mu / scale
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
dlaplace <- function(x,location=0,lambda=1,log=FALSE) {
    # f(y) = (1/(2b)) exp( -|y-a|/b )
    l <- (1/(2*lambda))*exp( -abs(x-location)/lambda )
    if ( log==TRUE ) l <- log(l)
    l
}

rlaplace <- function(n,location=0,lambda=1) {
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


dlkjcorr <- function( x , eta=1 , log=TRUE ) {
    det(x)^(eta-1)
}

rinvwishart <- function(s,df,Sigma,Prec) {
    #s-by-s Inverse Wishart matrix, df degree of freedom, precision matrix
    #Prec. Distribution of W^{-1} for Wishart W with nu=df+s-1 degree of
    # freedoom, covar martix Prec^{-1}.
    # NOTE  mean of riwish is proportional to Prec
    if ( missing(Prec) ) Prec <- solve(Sigma)
    if (df<=0) stop ("Inverse Wishart algorithm requires df>0")
    R <- diag(sqrt(2*rgamma(s,(df + s  - 1:s)/2)))
    R[outer(1:s, 1:s,  "<")] <- rnorm (s*(s-1)/2)
    S <- t(solve(R))%*% chol(Prec)
    return(t(S)%*%S)
}

log_sum_exp <- function( x ) {
    xmax <- max(x)
    xsum <- sum( exp( x - xmax ) )
    xmax + log(xsum)
}

# faster but less accurate version
#dzipois <- function( x , p , lambda , log=TRUE ) {
#    ll <- ifelse( x==0 , p + (1-p)*exp(-lambda) , (1-p)*dpois(x,lambda,FALSE) )
#    log(ll)
#}

# zero-inflated poisson distribution
# this is slow, but accurate
dzipois <- function(x,p,lambda,log=FALSE) {
    ll <- rep(0,length(x))
    p_i <- p[1]
    lambda_i <- lambda[1]
    for ( i in 1:length(x) ) {
        if ( length(p)>1 ) p_i <- p[i]
        if ( length(lambda)>1 ) lambda_i <- lambda[i]
        if ( x[i]==0 ) {
            ll[i] <- log_sum_exp( c( log(p_i) , log(1-p_i)+dpois(x[i],lambda_i,TRUE) ) )
        } else {
            ll[i] <- log(1-p_i) + dpois(x[i],lambda_i,TRUE)
        }
    }
    if ( log==FALSE ) ll <- exp(ll)
    return(ll)
}

rzipois <- function(n,p,lambda) {
    z <- rbinom(n,size=1,prob=p)
    y <- (1-z)*rpois(n,lambda)
    return(y)
}
