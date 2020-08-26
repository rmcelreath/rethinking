# distributions

# bernoulli convenience functions
dbern <- function(x,prob,log=FALSE) {
    dbinom(x,size=1,prob=prob,log=log)
}
rbern <- function(n,prob=0.5) {
    rbinom(n,size=1,prob=prob)
}

# generalized pareto
dgpareto <- function( x , u=0 , shape=1 , scale=1 , log=FALSE ) {
    if ( shape==0 ) {
        lp <- log(1/scale) - (x-u)/scale
    } else {
        lp <- log(1/scale) - (1/shape+1) * log(1 + shape*(x-u)/scale)
    }
    if ( log==FALSE ) lp <- exp(lp)
    return(lp)
}
rgpareto <- function( n , u=0 , shape=1 , scale=1 ) {
    x <- runif(n)
    if ( shape==0 )
        x <- u - scale * log(x)
    else
        x <- u + scale*(x^(-shape) - 1)/shape
    return(x)
}

dpareto <- function( x , xmin=0 , alpha=1 , log=FALSE ) {
    y <- dgpareto( x , u=xmin , shape=1/alpha , scale=xmin/alpha , log=log )
    return(y)
}
rpareto <- function( n , xmin=0 , alpha=1 ) {
    y <- rgpareto( n=n , u=xmin , shape=1/alpha , scale=xmin/alpha )
    return(y)
}

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

# complementary log log
# x is probability scale
cloglog <- function( x ) {
    log(-log(1-x))
}

# inverse cloglog
inv_cloglog <- function( x ) {
    p <- 1 - exp(-exp(x))
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

rgampois2 <- function (n , mu , scale) {
    p_p <- scale / (scale + mu)
    p_n <- scale
    rnbinom(n, size = p_n, prob = p_p)
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


# onion method correlation matrix - omits normalization constant
dlkjcorr <- function( x , eta=1 , log=TRUE ) {
    ll <- det(x)^(eta-1)
    if ( log==TRUE ) ll <- log(ll)
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
    if ( K==1 ) {
        X <- X[[1]]
        # vector or matrix?
        if ( is.null(dim(X)) ) {
            # vector
            f <- exp( X )
            denom <- sum(f)
            p <- as.numeric(f/denom)
            return(p)
        } else {
            # matrix - columns should be outcomes
            N <- nrow(X)
            f <- lapply( 1:N , function(i) exp(X[i,]) )
            denom <- sapply( 1:N , function(i) sum(f[[i]]) )
            p <- sapply( 1:N , function(i) unlist(f[[i]])/denom[i] )
            p <- t(as.matrix(p))
            return(p)
        }
    } else {
        # comma separated inputs
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
}


# dmvnorm2 - SRS form of multi_normal
# uses package mvtnorm
dmvnorm2 <- function( x , Mu , sigma , Rho , log=FALSE ) {
    DS <- diag(sigma)
    SIGMA <- DS %*% Rho %*% DS
    dmvnorm(x, Mu, SIGMA, log=log )
}

# vectorized version of rmvnorm, with split scale vector and correlation matrix
rmvnorm2 <- function( n , Mu=rep(0,length(sigma)) , sigma=rep(1,length(Mu)) , Rho=diag(length(Mu)) , method="chol" ) {
    ldim <- function(x) { 
        z <- length(dim(x))
        if ( z==0 ) z <- 1
        return(z)
    }
    if ( ldim(Mu)>1 || ldim(sigma)>1 || ldim(Rho)>2 ) {
        m <- 0 # dimension of multi_normal
        mR <- 0 # dim of Rho
        # vectorization needed
        # make Mu, sigma, and Rho same implied length
        K_Mu <- 1
        if ( ldim(Mu)==2 ) {
            K_Mu <- dim(Mu)[1] # matrix
            m <- dim(Mu)[2]
        } else {
            m <- length(Mu)
        }

        K_sigma <- 1
        if ( ldim(sigma)==2 ) {
            K_sigma <- dim(sigma)[1] # matrix
        }
        K_Rho <- 1
        if ( ldim(Rho)==3 ) {
            K_Rho <- dim(Rho)[1] # array (vector of matrices)
        }
        mR <- dim(Rho)[2] # should work whether Rho is matrix or array
        
        # get largest dim
        #K <- max( K_Mu , K_sigma , K_Rho )
        K <- n

        # need to fill to same size
        if ( K_Mu < K ) {
            Mu_new <- matrix( NA , nrow=K , ncol=m )
            row_list <- rep( 1:K_Mu , length.out=K )
            if ( K_Mu==1 )
                for ( i in 1:K )
                    Mu_new[i,] <- Mu
            else
                for ( i in 1:K )
                    Mu_new[i,] <- Mu[row_list[i],]
            Mu <- Mu_new
        }
        if ( K_sigma < K ) {
            sigma_new <- matrix( NA , nrow=K , ncol=m )
            row_list <- rep( 1:K_sigma , length.out=K )
            for ( i in 1:K )
                sigma_new[i,] <- sigma[row_list[i],]
            sigma <- sigma_new
        }
        if ( K_Rho < K ) {
            Rho_new <- array( NA , dim=c( K , m , m ) )
            row_list <- rep( 1:K_Rho , length.out=K )
            for ( i in 1:K )
                Rho_new[i,,] <- Rho[row_list[i],,]
            Rho <- Rho_new
        }
        # should all be same length now, so can brute force vectorize over and sample
        
        result <- array( NA , dim=c( K , m ) )
        for ( i in 1:n ) {
            DS <- diag(sigma[i,])
            SIGMA <- DS %*% Rho[i,,] %*% DS
            result[i,] <- rmvnorm( n=1 , mean=Mu[i,] , sigma=SIGMA , method=method )
        }
        return(result)

    } else {
        DS <- diag(sigma)
        SIGMA <- DS %*% Rho %*% DS
        return( rmvnorm( n=n , mean=Mu , sigma=SIGMA , method=method ) )
    }
}

# student t
dstudent <- function(x,nu=2,mu=0,sigma=1,log=FALSE) {
    #y <- gamma((nu+1)/2)/(gamma(nu/2)*sqrt(pi*nu)*sigma) * ( 1 + (1/nu)*((x-mu)/sigma)^2 )^(-(nu+1)/2)
    y <- lgamma((nu+1)/2) - lgamma(nu/2) - 0.5*log(pi*nu) - log(sigma) + (-(nu+1)/2)*log( 1 + (1/nu)*((x-mu)/sigma)^2 )
    if ( log==FALSE ) y <- exp(y)
    return(y)
}

rstudent <- function(n,nu=2,mu=0,sigma=1) {
    y <- rt(n,df=nu)
    y <- y*sigma + mu
    return(y)
}
