# distributions

# ordered categorical density functions
#' @export
logistic <- function( x ) {
    p <- 1 / ( 1 + exp( -x ) )
    p <- ifelse( x==Inf , 1 , p )
    p
}

#' @export
inv_logit <- function( x ) {
    p <- 1 / ( 1 + exp( -x ) )
    p <- ifelse( x==Inf , 1 , p )
    p
}

#' @export
logit <- function( x ) {
    log( x ) - log( 1-x )
}

#' @export
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



#' Ordered categorical log-odds probability density
#'
#' Functions for computing density, cumulative probability and producing random
#' samples from an ordered categorical probability density.
#'
#' These functions provide for common density calculations for the ordered
#' categorical probability density commonly known as an "ordered logit" or
#' "ordered logistic." This is the same probability density assumed by the
#' \code{polr} function (library \code{MASS}).
#'
#' @aliases dordlogit rordlogit pordlogit
#' @param x Integer values to compute densities or probabilies of
#' @param a Vector of log-odds intercepts
#' @param phi Linear model of log-odds
#' @param log If \code{TRUE}, returns log-probability instead of probability
#' @author Richard McElreath
#' @export
dordlogit <- function( x , phi , a , log=FALSE ) {
    a <- c( as.numeric(a) , Inf )
    p <- logistic( a[x] - phi )
    na <- c( -Inf , a )
    np <- logistic( na[x] - phi )
    p <- p - np
    if ( log==TRUE ) p <- log(p)
    p
}

#' @export
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

#' @export
dbeta2 <- function( x , prob , theta , log=FALSE ) {
    a <- prob * theta
    b <- (1-prob) * theta
    dbeta( x , shape1=a , shape2=b , log=log )
}

#' @export
rbeta2 <- function( n , prob , theta ) {
    a <- prob * theta
    b <- (1-prob) * theta
    rbeta( n , shape1=a , shape2=b )
}

# Ben Bolker's dbetabinom from emdbook package


#' Beta-binomial probability density
#'
#' Functions for computing density and producing random samples from a
#' beta-binomial probability distribution.
#'
#' These functions provide density and random number calculations for
#' beta-binomial observations. The \code{dbetabinom} code is based on Ben
#' Bolker's original code in the \code{emdbook} package.
#'
#' Either \code{prob} and \code{theta} OR \code{shape1} and \code{shape2} must
#' be provided. The two parameterizations are related by shape1 = prob * theta,
#' shape2 = (1-prob) * theta.
#'
#' The \code{rbetabinom} function generates random beta-binomial observations
#' by using both \code{\link{rbeta}} and \code{\link{rbinom}}. It draws
#' \code{n} values from a beta distribution. Then for each, it generates a
#' random binomial observation.
#'
#' @aliases dbetabinom rbetabinom
#' @param x Integer values to compute probabilies of
#' @param size Number of trials
#' @param prob Average probability of beta distribution
#' @param theta Dispersion of beta distribution
#' @param shape1 First shape parameter of beta distribution (alpha)
#' @param shape2 Second shape parameter of beta distribution (beta)
#' @param log If \code{TRUE}, returns log-probability instead of probability
#' @param n Number of random observations to sample
#' @author Richard McElreath
#' @export
#' @examples
#'
#' \dontrun{
#' data(reedfrogs)
#' reedfrogs$pred_yes <- ifelse( reedfrogs$pred=="pred" , 1 , 0 )
#'
#' # map model fit
#' # note exp(log_theta) to constrain theta to positive reals
#' m <- map(
#'     alist(
#'         surv ~ dbetabinom( density , p , exp(log_theta) ),
#'         logit(p) <- a + b*pred_yes,
#'         a ~ dnorm(0,10),
#'         b ~ dnorm(0,1),
#'         log_theta ~ dnorm(1,10)
#'     ),
#'     data=reedfrogs )
#'
#' # map2stan model fit
#' # constraint on theta is passed via contraints list
#' m.stan <- map2stan(
#'     alist(
#'         surv ~ dbetabinom( density , p , theta ),
#'         logit(p) <- a + b*pred_yes,
#'         a ~ dnorm(0,10),
#'         b ~ dnorm(0,1),
#'         theta ~ dcauchy(0,1)
#'     ),
#'     data=reedfrogs,
#'     constraints=list(theta="lower=0") )
#' }
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

#' @export
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
#' @export
dgamma2 <- function( x , mu , scale , log=FALSE ) {
    dgamma( x , shape=mu / scale , scale=scale , log=log )
}

#' @export
rgamma2 <- function( n , mu , scale ) {
    rgamma( n , shape=mu/scale , scale=scale )
}



#' Gamma-Poisson probability density
#'
#' Functions for computing density and producing random samples from a
#' gamma-Poisson (negative-binomial) probability distribution.
#'
#' These functions provide density and random number calculations for
#' gamma-Poisson observations. These functions use \code{dnbinom} and
#' \code{rnbinom} internally, but convert the parameters from the \code{mu} and
#' \code{scale} form. The \code{size} parameter of the negative-binomial is
#' defined by \code{mu/scale}, and the \code{prob} parameter of the
#' negative-binomial is the same as \code{1/(1+scale)}.
#'
#' @aliases dgampois rgampois
#' @param x Integer values to compute probabilies of
#' @param mu Mean of gamma distribution
#' @param scale Scale parameter of gamma distribution
#' @param log If \code{TRUE}, returns log-probability instead of probability
#' @param n Number of random observations to sample
#' @author Richard McElreath
#' @export
#' @examples
#'
#' \dontrun{
#' data(Hurricanes)
#'
#' # map model fit
#' # note exp(log_scale) to constrain scale to positive reals
#' m <- map(
#'     alist(
#'         deaths ~ dgampois( mu , exp(log_scale) ),
#'         log(mu) <- a + b*femininity,
#'         a ~ dnorm(0,100),
#'         b ~ dnorm(0,1),
#'         log_scale ~ dnorm(1,10)
#'     ),
#'     data=Hurricanes )
#'
#' # map2stan model fit
#' # constraint on scale is passed via contraints list
#' m.stan <- map2stan(
#'     alist(
#'         deaths ~ dgampois( mu , scale ),
#'         log(mu) <- a + b*femininity,
#'         a ~ dnorm(0,100),
#'         b ~ dnorm(0,1),
#'         scale ~ dcauchy(0,2)
#'     ),
#'     data=Hurricanes,
#'     constraints=list(scale="lower=0"),
#'     start=list(scale=2) )
#' }
dgampois <- function( x , mu , scale , log=FALSE ) {
    shape <- mu / scale
    prob <- 1 / ( 1 + scale )
    dnbinom( x , size=shape , prob=prob , log=log )
}

#' @export
rgampois <- function( n , mu , scale ) {
    shape <- mu / scale
    prob <- 1 / ( 1 + scale )
    rnbinom( n , size=shape , prob=prob )
}

# laplace (double exponential)
#' @export
dlaplace <- function(x,location=0,lambda=1,log=FALSE) {
    # f(y) = (1/(2b)) exp( -|y-a|/b )
    # l <- (1/(2*lambda))*exp( -abs(x-location)/lambda )
    ll <- -abs(x-location)/lambda - log(2*lambda)
    if ( log==FALSE ) ll <- exp(ll)
    ll
}

#' @export
rlaplace <- function(n,location=0,lambda=1) {
    # generate from difference of two random exponentials with rate 1/lambda
    y <- rexp(n=n,rate=1/lambda) - rexp(n=n,rate=1/lambda)
    y <- y + location # offset location
    y
}


# onion method correlation matrix


#' LKJ correlation matrix probability density
#'
#' Functions for computing density and producing random samples from the LKJ
#' onion method correlation matrix distribution.
#'
#' The LKJ correlation matrix distribution is based upon work by Lewandowski,
#' Kurowicka, and Joe. When the parameter \code{eta} is equal to 1, it defines
#' a flat distribution of correlation matrices. When \code{eta > 1}, the
#' distribution is instead concentrated towards to identity matrix. When
#' \code{eta < 1}, the distribution is more concentrated towards extreme
#' correlations at -1 or +1.
#'
#' It can be easier to understand this distribution if we recognize that the
#' individual correlations within the matrix follow a beta distribution defined
#' on -1 to +1. Thus \code{eta} resembles \code{theta} in the beta
#' parameterization with a mean p and scale (sample size) theta.
#'
#' The function \code{rlkjcorr} returns an 3-dimensional array in which the
#' first dimension indexes matrices. In the event that \code{n=1}, it returns
#' instead a single matrix.
#'
#' @aliases dlkjcorr rlkjcorr
#' @param x Matrix to compute probability density for
#' @param eta Parameter controlling shape of distribution
#' @param K Dimension of correlation matrix
#' @param log If \code{TRUE}, returns log-probability instead of probability
#' @param n Number of random matrices to sample
#' @author Richard McElreath
#' @references Lewandowski, Kurowicka, and Joe. 2009. Generating random
#' correlation matrices based on vines and extended onion method. Journal of
#' Multivariate Analysis. 100:1989-2001.
#'
#' Stan Modeling Language User's Guide and Reference Manual, v2.6.2
#'
#' @export
#' @examples
#' R <- rlkjcorr(n=1,K=2,eta=4)
#' dlkjcorr(R,4)
#'
#' # plot density of correlation
#' R <- rlkjcorr(1e4,K=2,eta=4)
#' dens( R[,1,2] )
#'
#' # visualize 3x3 matrix
#' R <- rlkjcorr(1e3,K=3,eta=2)
#' plot( R[,1,2] , R[,1,3] , col=col.alpha("black",0.2) , pch=16 )
#'
dlkjcorr <- function( x , eta=1 , log=TRUE ) {
    ll <- det(x)^(eta-1)
    if ( log==FALSE ) ll <- exp(ll)
    return(ll)
}

# Ben's rLKJ function
# K is dimension of matrix
#' @export
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
#' @export
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



#' Rethinking utility functions
#'
#' Numerically stable computation of log sums of exponentiated values.
#'
#' This function is used by \code{WAIC} to compute the log average probability
#' used in the formula for WAIC. In that context, a vector of log-probabilities
#' is passed to \code{log_sum_exp}, obviating the need to explicitly
#' exponentiate. This helps to avoid rounding errors that occur when working
#' with direct probabilities.
#'
#' @param x vector of values
#' @author Richard McElreath
#' @examples
#'
#'
#' @export
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


#' Zero-inflated Poisson probability density
#'
#' Functions for computing density and producing random samples from a
#' zero-inflated Poisson probability distribution.
#'
#' These functions provide density and random number calculations for
#' zero-inflated Poisson observations. This distribution is defined as a finite
#' mixture of zeros and Poisson values, where \code{p} determines the weight of
#' the pure zeros. As such, the probability of a zero is \code{p +
#' (1-p)exp(-lambda)}. And the probability of a non-zero value \code{x} is
#' \code{(1-p)lambda^x exp(-lambda)/x!}.
#'
#' \code{dzipois} does its calculations in a way meant to reduce numerical
#' issues with probabilities close to 0 and 1. As a result, it's not very fast.
#'
#' @aliases dzipois rzipois
#' @param x Integer values to compute densities or probabilies of
#' @param p Probability of a zero from bernoulli process
#' @param lambda Poisson rate parameter (mean)
#' @param log If \code{TRUE}, returns log-probability instead of probability
#' @author Richard McElreath
#' @export
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

#' @export
rzipois <- function(n,p,lambda) {
    z <- rbinom(n,size=1,prob=p)
    y <- (1-z)*rpois(n,lambda)
    return(y)
}

# zero-augmented gamma2 distribution


#' Zero-augmented gamma probability density
#'
#' Function for computing density from a zero-augmented gamma probability
#' distribution.
#'
#' This distribution is defined as a finite mixture of zeros and strictly
#' positive gamma distributed values, where \code{prob} determines the weight
#' of the zeros. As such, the probability of a zero is \code{prob}, and the
#' probability of a non-zero value \code{x} is \code{(1-prob)*dgamma( x ,
#' mu/scale , scale )}.
#'
#' @param x Values to compute densities for
#' @param prob Probability of a zero
#' @param mu Mean of gamma distribution
#' @param scale Scale parameter (same as \code{scale} in \code{\link{dgamma}})
#' @param log If \code{TRUE}, returns log-density instead of density
#' @author Richard McElreath
#' @export
#' @examples
#'
dzagamma2 <- function(x,prob,mu,scale,log=FALSE) {
    K <- as.data.frame(cbind(x=x,prob=prob,mu=mu,scale=scale))
    llg <- dgamma( x , shape=mu/scale , scale=scale , log=TRUE )
    ll <- ifelse( K$x==0 , log(K$prob) , log(1-K$prob) + llg )
    if ( log==FALSE ) ll <- exp(ll)
    ll
}

# zero-inflated binomial distribution
# this is slow, but accurate


#' Zero-inflated binomial probability density
#'
#' Functions for computing density and producing random samples from a
#' zero-inflated binomial probability distribution.
#'
#' These functions provide density and random number calculations for
#' zero-inflated binomial observations. This distribution is defined as a
#' finite mixture of zeros and binomial values, where \code{p_zero} determines
#' the weight of the pure zeros. As such, the probability of a zero is
#' \code{p_zero + (1-p_zero)(1-prob)^size}. And the probability of a non-zero
#' value \code{x} is \code{(1-p_zero) choose(size,x) prob^x (1-prob)^(size-x)}.
#'
#' \code{dzibinom} does its calculations in a way meant to reduce numerical
#' issues with probabilities close to 0 and 1. As a result, it's not very fast.
#'
#' @aliases dzibinom rzibinom
#' @param x Integer values to compute densities or probabilies of
#' @param p_zero Probability of a zero from bernoulli process
#' @param size Number of binomial trials
#' @param prob Probability of success in binomial trial
#' @param log If \code{TRUE}, returns log-probability instead of probability
#' @author Richard McElreath
#' @export
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

#' @export
rzibinom <- function(n,p_zero,size,prob) {
    z <- rbinom(n,size=1,prob=p_zero)
    y <- (1-z)*rbinom(n,size,prob)
    return(y)
}

# Student's t, with scale and location parameters

#dstudent <- function( x , nu , mu , sigma , log=TRUE ) {
#
#}

#' @export
# generalized normal
dgnorm <- function( x , mu , alpha , beta , log=FALSE ) {
    ll <- log(beta) - log( 2*alpha*gamma(1/beta) ) - (abs(x-mu)/alpha)^beta
    if ( log==FALSE ) ll <- exp(ll)
    return(ll)
}

# categorical distribution for multinomial models

#' @export
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

#' @export
rcategorical <- function( n , prob ) {
    k <- length(prob)
    y <- sample( 1:k , size=n , prob=prob , replace=TRUE )
    return(y)
}

### softmax rule for multinomial
#' @export
multilogistic <- function( x , lambda=1 , diff=TRUE , log=FALSE ) {
    if ( diff==TRUE ) x <- x - min(x)
    f <- exp( lambda * x )
    if ( log==FALSE ) result <- f / sum(f)
    if ( log==TRUE ) result <- log(f) - log(sum(f))
    result
}

# multinomial logit link
# needs to recognize vectors, as well
#' @export
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


#' Multivariate Gaussian probability density
#'
#' This is an alternative parameterization of the ordinary multivariate
#' Gaussian probability density.
#'
#' These functions merely compose the variance-covariance matrix from separate
#' standard deviation and correlation matrix arguments. They then use
#' \code{dmvnorm} and \code{rmvnorm} from the \code{mvtnorm} package to perform
#' calculations.
#'
#' @aliases dmvnorm2 rmvnorm2
#' @param x Values to compute densities of
#' @param Mu Mean vector
#' @param sigma Vector of standard deviations
#' @param Rho Correlation matrix
#' @param log If \code{TRUE}, returns log-density instead of density
#' @param n Number of random observations to sample
#' @author Richard McElreath
#' @export
#' @examples
#'
#' dmvnorm2( c(1,0) , Mu=c(0,0) , sigma=c(1,1) , Rho=diag(2) )
#' rmvnorm2( 10 , Mu=c(1,2) )
#'
dmvnorm2 <- function( x , Mu , sigma , Rho , log=FALSE ) {
    DS <- diag(sigma)
    SIGMA <- DS %*% Rho %*% DS
    dmvnorm(x, Mu, SIGMA, log=log )
}

#' @export
rmvnorm2 <- function( n , Mu=rep(0,length(sigma)) , sigma=rep(1,length(Mu)) , Rho=diag(length(Mu)) , method="chol" ) {
    DS <- diag(sigma)
    SIGMA <- DS %*% Rho %*% DS
    rmvnorm( n=n , mean=Mu , sigma=SIGMA , method=method )
}

