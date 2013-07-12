# compare

# compare class definition and show method
setClass( "compareIC" , representation( output="data.frame" ) )
compare.show <- function( object ) {
    print( round( object@output , 2 ) )
}
setMethod( "show" , "compareIC" , function(object) compare.show(object) )

# AICc/BIC model comparison table
compare <- function( ... , nobs=NULL , sort="AICc" , delta=TRUE , digits=4 ) {
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
    BIC.list <- sapply( L , function(z) myBIC( z , nobs=nobs ) )
    dAICc <- AICc.list - min( AICc.list )
    dBIC <- BIC.list - min( BIC.list )
    
    post.AICc <- exp( -0.5*dAICc ) / sum( exp(-0.5*dAICc) )
    post.BIC <- exp( -0.5*dBIC ) / sum( exp(-0.5*dBIC) )
    
    #post.AICc <- format.pval( post.AICc , digits=digits )
    #post.BIC <- format.pval( post.BIC , digits=digits )
    
    k <- sapply( L , getdf )
    
    result <- data.frame( k=k , AICc=AICc.list , BIC=BIC.list , w.AICc=post.AICc , w.BIC=post.BIC )
    if ( delta==TRUE ) {
        r2 <- data.frame( dAICc=dAICc , dBIC=dBIC )
        result <- cbind( result , r2 )
    }
    rownames( result ) <- mnames
    
    if ( !is.null(sort) ) {
        result <- result[ order( result$AICc ) , ]
    }
    
    new( "compareIC" , output=result )
}

