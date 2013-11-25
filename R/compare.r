# compare

# compare class definition and show method
setClass( "compareIC" , representation( output="data.frame" ) )
compare.show <- function( object ) {
    print( round( object@output , 2 ) )
}
setMethod( "show" , "compareIC" , function(object) compare.show(object) )

# AICc/BIC model comparison table
compare <- function( ... , nobs=NULL , sort="AICc" , BIC=FALSE , DIC=FALSE , delta=TRUE , DICsamples=1e4 ) {
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
                dev <- sapply( 1:nrow(post) , function(i) 2*m@fminuslogl( post[i,] ) )
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

