# distributions

# ordered categorical density functions
logistic <- function( x ) {
    p <- 1 / ( 1 + exp( -x ) )
    p <- ifelse( x==Inf , 1 , p )
    p
}

inv_logit <- function( x ) {
    p <- 1 / ( 1 + exp( -x ) )
    p <- ifelse( x==Inf , 1 , p )
    p
}

logit <- function( x ) {
    log( x ) - log( 1-x )
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
    if ( length(phi)==1 ) {
        p <- dordlogit( k , a=a , phi=phi , log=FALSE )
        y <- sample( k , size=n , replace=TRUE , prob=p )
    } else {
        # vectorized input
        y <- rep(NA,n)
        if ( n > length(phi) ) {
            # need to repeat some phi values
            phi <- rep(phi, ceiling(n/length(phi)) )
        }
        for ( i in 1:n ) {
            p <- dordlogit( k , a=a , phi=phi[i] , log=FALSE )
            y[i] <- sample( k , size=1 , replace=TRUE , prob=p )
        }
    }
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

rbetabinom <- function( n , size, prob, theta, shape1, shape2 ) {
    if (missing(prob) && !missing(shape1) && !missing(shape2)) {
        prob <- shape1/(shape1 + shape2)
        theta <- shape1 + shape2
    }
    # first generate beta probs
    p <- rbeta2( n , prob , theta )
    # then generate binomial counts
    y <- rbinom( n , size , p )
    return(y)
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

# laplace (double exponential)
dlaplace <- function(x,location=0,lambda=1,log=FALSE) {
    # f(y) = (1/(2b)) exp( -|y-a|/b )
    # l <- (1/(2*lambda))*exp( -abs(x-location)/lambda )
    ll <- -abs(x-location)/lambda - log(2*lambda)
    if ( log==FALSE ) ll <- exp(ll)
    ll
}

rlaplace <- function(n,location=0,lambda=1) {
    # generate from difference of two random exponentials with rate 1/lambda
    y <- rexp(n=n,rate=1/lambda) - rexp(n=n,rate=1/lambda)
    y <- y + location # offset location
    y
}


# onion method correlation matrix
dlkjcorr <- function( x , eta=1 , log=TRUE ) {
    ll <- det(x)^(eta-1)
    if ( log==FALSE ) ll <- exp(ll)
    return(ll)
}

# Ben's rLKJ function
# K is dimension of matrix
rlkjcorr <- function ( n , K , eta = 1 ) {

    stopifnot(is.numeric(K), K >= 2, K == as.integer(K))
    stopifnot(eta > 0)
    #if (K == 1) return(matrix(1, 1, 1))
    
    f <- function() {
        alpha <- eta + (K - 2)/2
        r12 <- 2 * rbeta(1, alpha, alpha) - 1
        R <- matrix(0, K, K) # upper triangular Cholesky factor until return()
        R[1,1] <- 1
        R[1,2] <- r12
        R[2,2] <- sqrt(1 - r12^2)
        if(K > 2) for (m in 2:(K - 1)) {
            alpha <- alpha - 0.5
            y <- rbeta(1, m / 2, alpha)
            
            # Draw uniformally on a hypersphere
            z <- rnorm(m, 0, 1)
            z <- z / sqrt(crossprod(z)[1])
            
            R[1:m,m+1] <- sqrt(y) * z
            R[m+1,m+1] <- sqrt(1 - y)
        }
        return(crossprod(R))
    }
    R <- replicate( n , f() )
    if ( dim(R)[3]==1 ) {
        R <- R[,,1]
    } else {
        # need to move 3rd dimension to front, so conforms to array structure that Stan uses
        R <- aperm(R,c(3,1,2))
    }
    return(R)
}

# wishart
rinvwishart <- function(s,df,Sigma,Prec) {
    #s-by-s Inverse Wishart matrix, df degree of freedom, precision matrix
    #Prec. Distribution of W^{-1} for Wishart W with nu=df+s-1 degree of
    # freedom, covar martix Prec^{-1}.
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

# zero-augmented gamma2 distribution
dzagamma2 <- function(x,prob,mu,scale,log=FALSE) {
    K <- as.data.frame(cbind(x=x,prob=prob,mu=mu,scale=scale))
    llg <- dgamma( x , shape=mu/scale , scale=scale , log=TRUE )
    ll <- ifelse( K$x==0 , log(K$prob) , log(1-K$prob) + llg )
    if ( log==FALSE ) ll <- exp(ll)
    ll
}

# zero-inflated binomial distribution
# this is slow, but accurate
dzibinom <- function(x,p_zero,size,prob,log=FALSE) {
    ll <- rep(0,length(x))
    pz_i <- p_zero[1]
    size_i <- size[1]
    prob_i <- prob[1]
    for ( i in 1:length(x) ) {
        if ( length(p_zero)>1 ) pz_i <- p_zero[i]
        if ( length(size)>1 ) size_i <- size[i]
        if ( length(prob)>1 ) prob_i <- prob[i]
        if ( x[i]==0 ) {
            ll[i] <- log_sum_exp( c( log(pz_i) , log(1-pz_i)+dbinom(x[i],size_i,prob_i,TRUE) ) )
        } else {
            ll[i] <- log(1-pz_i) + dbinom(x[i],size_i,prob_i,TRUE)
        }
    }
    if ( log==FALSE ) ll <- exp(ll)
    return(ll)
}

rzibinom <- function(n,p_zero,size,prob) {
    z <- rbinom(n,size=1,prob=p_zero)
    y <- (1-z)*rbinom(n,size,prob)
    return(y)
}

# Student's t, with scale and location parameters

#dstudent <- function( x , nu , mu , sigma , log=TRUE ) {
#    
#}

# generalized normal
dgnorm <- function( x , mu , alpha , beta , log=FALSE ) {
    ll <- log(beta) - log( 2*alpha*gamma(1/beta) ) - (abs(x-mu)/alpha)^beta
    if ( log==FALSE ) ll <- exp(ll)
    return(ll)
}

# categorical distribution for multinomial models

dcategorical <- function( x , prob , log=TRUE ) {
    if ( class(prob)=="matrix" ) {
        # vectorized probability matrix
        # length of x needs to match nrow(prob)
        logp <- sapply( 1:nrow(prob) , function(i) log(prob[i,x[i]]) )
        if ( log==FALSE ) logp <- exp(logp)
        return( logp )
    } else {
        logp <- log(prob[x])
        if ( log==FALSE ) logp <- exp(logp)
        return(logp)
    }
}

rcategorical <- function( n , prob ) {
    k <- length(prob)
    y <- sample( 1:k , size=n , prob=prob , replace=TRUE )
    return(y)
}

### softmax rule for multinomial
multilogistic <- function( x , lambda=1 , diff=TRUE , log=FALSE ) {
    if ( diff==TRUE ) x <- x - min(x)
    f <- exp( lambda * x )
    if ( log==FALSE ) result <- f / sum(f)
    if ( log==TRUE ) result <- log(f) - log(sum(f))
    result
}

# multinomial logit link
# needs to recognize vectors, as well
softmax <- function( ... ) {
    X <- list(...)
    K <- length(X)
    X <- as.data.frame(X)
    N <- nrow(X)
    if ( N==1 ) {
        # just a vector of probabilities
        f <- exp( X[1,] )
        denom <- sum(f)
        p <- as.numeric(f/denom)
        return(p)
    } else {
        f <- lapply( 1:N , function(i) exp(X[i,]) )
        denom <- sapply( 1:N , function(i) sum(f[[i]]) )
        p <- sapply( 1:N , function(i) unlist(f[[i]])/denom[i] )
        p <- t(as.matrix(p))
        colnames(p) <- NULL
        return(p)
    }
}


# dmvnorm2 - SRS form of multi_normal
# uses package mvtnorm
dmvnorm2 <- function( x , Mu , sigma , Rho , log=FALSE ) {
    DS <- diag(sigma)
    SIGMA <- DS %*% Rho %*% DS
    dmvnorm(x, Mu, SIGMA, log=log )
}

rmvnorm2 <- function( n , Mu=rep(0,length(sigma)) , sigma=rep(1,length(Mu)) , Rho=diag(length(Mu)) , method="chol" ) {
    DS <- diag(sigma)
    SIGMA <- DS %*% Rho %*% DS
    rmvnorm( n=n , mean=Mu , sigma=SIGMA , method=method )
}

