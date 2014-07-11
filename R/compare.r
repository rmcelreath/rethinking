# compare

# compare class definition and show method
setClass( "compareIC" , representation( output="data.frame" ) )
compare.show <- function( object ) {
    print( round( object@output , 2 ) )
}
setMethod( "show" , "compareIC" , function(object) compare.show(object) )

# new compare function, defaulting to DIC
compare <- function( ... , n=1e3 , sort="DIC" , WAIC=FALSE ) {
    # retrieve list of models
    L <- list(...)
    if ( is.list(L[[1]]) && length(L)==1 )
        L <- L[[1]]
    
    # retrieve model names from function call
    mnames <- match.call()
    mnames <- as.character(mnames)[2:(length(L)+1)]
    
    if ( WAIC==FALSE ) {
        DIC.list <- lapply( L , function(z) DIC( z , n=n ) )
        pD.list <- sapply( DIC.list , function(x) attr(x,"pD") )
    } else {
        # use WAIC instead of DIC
        WAIC.list <- lapply( L , function(z) WAIC( z , n=n ) )
        pD.list <- sapply( WAIC.list , function(x) attr(x,"pWAIC") )
        DIC.list <- WAIC.list
    }
    
    DIC.list <- unlist(DIC.list)
    
    dDIC <- DIC.list - min( DIC.list )
    w.DIC <- ICweights( DIC.list )
    
    if ( WAIC==FALSE )
        result <- data.frame( DIC=DIC.list , pD=pD.list , dDIC=dDIC , weight=w.DIC )
    else
        result <- data.frame( WAIC=DIC.list , pWAIC=pD.list , dWAIC=dDIC , weight=w.DIC )
    
    rownames(result) <- mnames
    
    if ( !is.null(sort) ) {
        if ( WAIC==TRUE & sort=="DIC" ) sort <- "WAIC"
        result <- result[ order( result[[sort]] ) , ]
    }
    
    new( "compareIC" , output=result )
}

# plot method for compareIC results shows deviance in and expected deviance out of sample, for each model, ordered top-to-bottom by rank
setMethod("plot" , "compareIC" , function(x,y,...) {
    dev_in <- x@output[[1]] - x@output[[2]]
    dev_out <- x@output[[1]]
    n <- length(dev_in)
    dotchart( dev_in[n:1] , labels=rownames(x@output)[n:1] , xlab="deviance" , pch=16 , xlim=c(min(dev_in),max(dev_out)) , ... )
    points( dev_out[n:1] , 1:n )
})

# AICc/BIC model comparison table
compare_old <- function( ... , nobs=NULL , sort="AICc" , BIC=FALSE , DIC=FALSE , delta=TRUE , DICsamples=1e4 ) {
    require(bbmle)
    
    if ( is.null(nobs) ) {
        stop( "Must specify number of observations (nobs)." )
    }
    
    getdf <- function(x) {
        if (!is.null(df <- attr(x, "df"))) 
            return(df)
        else if (!is.null(df <- attr(logLik(x), "df"))) 
            return(df)
    }
    
    # need own BIC, as one in stats doesn't allow nobs
    myBIC <- function(x,nobs) {
        k <- getdf(x)
        as.numeric( -2*logLik(x) + log(nobs)*k )
    }
    
    # retrieve list of models
    L <- list(...)
    if ( is.list(L[[1]]) && length(L)==1 )
        L <- L[[1]]
    
    # retrieve model names from function call
    mnames <- match.call()
    mnames <- as.character(mnames)[2:(length(L)+1)]
    
    AICc.list <- sapply( L , function(z) AICc( z , nobs=nobs ) )
    dAICc <- AICc.list - min( AICc.list )
    post.AICc <- exp( -0.5*dAICc ) / sum( exp(-0.5*dAICc) )
    if ( BIC==TRUE ) {
        BIC.list <- sapply( L , function(z) myBIC( z , nobs=nobs ) )
        dBIC <- BIC.list - min( BIC.list )
        post.BIC <- exp( -0.5*dBIC ) / sum( exp(-0.5*dBIC) )
    }
    
    k <- sapply( L , getdf )
    
    result <- data.frame( k=k , AICc=AICc.list , w.AICc=post.AICc )
    if ( BIC==TRUE ) 
        result <- data.frame( k=k , AICc=AICc.list , BIC=BIC.list , w.AICc=post.AICc , w.BIC=post.BIC )

    if ( delta==TRUE ) {
        r2 <- data.frame( dAICc=dAICc )
        if ( BIC==TRUE ) r2 <- data.frame( dAICc=dAICc , dBIC=dBIC )
        result <- cbind( result , r2 )
    }

    # DIC from quadratic approx posterior defined by vcov and coef
    if ( DIC==TRUE ) {
        DIC.list <- rep( NA , length(L) )
        pD.list <- rep( NA , length(L) )
        for ( i in 1:length(L) ) {
            m <- L[[i]]
            if ( class(m)=="map" ) {
                post <- sample.qa.posterior( m , n=DICsamples )
                message( paste("Computing DIC for model",mnames[i]) )
                dev <- sapply( 1:nrow(post) , 
                    function(i) {
                        p <- post[i,]
                        names(p) <- names(post)
                        2*m@fminuslogl( p ) 
                    }
                )
                dev.hat <- deviance(m)
                DIC.list[i] <- dev.hat + 2*( mean(dev) - dev.hat )
                pD.list[i] <- ( DIC.list[i] - dev.hat )/2
            }
        }
        ddic <- DIC.list - min(DIC.list)
        wdic <- exp( -0.5*ddic ) / sum( exp(-0.5*ddic) )
        rdic <- data.frame( DIC=as.numeric(DIC.list) , pD=pD.list , wDIC=wdic , dDIC=ddic )
        result <- cbind( result , rdic )
    }

    # add model names to rows
    rownames( result ) <- mnames
    
    if ( !is.null(sort) ) {
        result <- result[ order( result[[sort]] ) , ]
    }
    
    new( "compareIC" , output=result )
}

# convert estimated D_test values to weights
ICweights <- function( dev ) {
    d <- dev - min(dev)
    f <- exp(-0.5*d)
    w <- f/sum(f)
    return(w)
}

# build ensemble of samples using DIC/WAIC weights
ensemble <- function( ... ) {
    L <- list(...)
}

# tests
if (FALSE) {

library(rethinking)
data(chimpanzees)

d <- list( 
    pulled_left = chimpanzees$pulled.left ,
    prosoc_left = chimpanzees$prosoc.left ,
    condition = chimpanzees$condition ,
    actor = as.integer( chimpanzees$actor )
)

m0 <- map(
    alist(
        pulled_left ~ dbinom(1,theta),
        logit(theta) <- a,
        a ~ dnorm(0,1)
    ) ,
    data=d,
    start=list(a=0)
)

m1 <- map(
    alist(
        pulled_left ~ dbinom(1,theta),
        logit(theta) <- a + bp*prosoc_left,
        a ~ dnorm(0,1),
        bp ~ dnorm(0,1)
    ) ,
    data=d,
    start=list(a=0,bp=0)
)

m2 <- map(
    alist(
        pulled_left ~ dbinom(1,theta),
        logit(theta) <- a + bp*prosoc_left + bpc*condition*prosoc_left,
        a ~ dnorm(0,1),
        bp ~ dnorm(0,1),
        bpc ~ dnorm(0,1)
    ) ,
    data=d,
    start=list(a=0,bp=0,bpc=0)
)

m3 <- map(
    alist(
        pulled_left ~ dbinom(1,theta),
        logit(theta) <- a + bp*prosoc_left + bc*condition + bpc*condition*prosoc_left,
        a ~ dnorm(0,1),
        bp ~ dnorm(0,1),
        bc ~ dnorm(0,1),
        bpc ~ dnorm(0,1)
    ) ,
    data=d,
    start=list(a=0,bp=0,bc=0,bpc=0)
)

( x <- compare(m0,m1,m2,m3) )

plot(x)

# now map2stan

m0 <- map2stan(
    alist(
        pulled_left ~ dbinom(1,theta),
        logit(theta) <- a,
        a ~ dnorm(0,1)
    ) ,
    data=d,
    start=list(a=0)
)

m1 <- map2stan(
    alist(
        pulled_left ~ dbinom(1,theta),
        logit(theta) <- a + bp*prosoc_left,
        a ~ dnorm(0,1),
        bp ~ dnorm(0,1)
    ) ,
    data=d,
    start=list(a=0,bp=0)
)

m2 <- map2stan(
    alist(
        pulled_left ~ dbinom(1,theta),
        logit(theta) <- a + bp*prosoc_left + bpc*condition*prosoc_left,
        a ~ dnorm(0,1),
        bp ~ dnorm(0,1),
        bpc ~ dnorm(0,1)
    ) ,
    data=d,
    start=list(a=0,bp=0,bpc=0)
)

m3 <- map2stan(
    alist(
        pulled_left ~ dbinom(1,theta),
        logit(theta) <- a + bp*prosoc_left + bc*condition + bpc*condition*prosoc_left,
        a ~ dnorm(0,1),
        bp ~ dnorm(0,1),
        bc ~ dnorm(0,1),
        bpc ~ dnorm(0,1)
    ) ,
    data=d,
    start=list(a=0,bp=0,bc=0,bpc=0)
)

m4 <- map2stan(
    alist(
        pulled_left ~ dbinom(1,theta),
        logit(theta) <- a + aj + bp*prosoc_left + bc*condition + bpc*condition*prosoc_left,
        a ~ dnorm(0,1),
        aj[actor] ~ dnorm(0,sigma_actor),
        bp ~ dnorm(0,1),
        bc ~ dnorm(0,1),
        bpc ~ dnorm(0,1),
        sigma_actor ~ dcauchy(0,1)
    ) ,
    data=d,
    start=list(a=0,bp=0,bc=0,bpc=0,sigma_actor=1,aj=rep(0,7))
)

( x1 <- compare(m0,m1,m2,m3,m4,WAIC=FALSE) )

( x2 <- compare(m0,m1,m2,m3,m4,WAIC=TRUE) )

plot(x1)

plot(x2)


}
