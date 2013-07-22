# coeftab

# coeftab class definition and show method
setClass( "coeftab" , representation( coefs="matrix" , nobs="numeric" , AIC="numeric" , digits="numeric" , width="numeric" ) )
coeftab.show <- function( object ) {
    result <- object@coefs
    if ( !is.null(object@nobs) ) {
        result <- rbind( result , object@nobs )
        rownames(result)[ nrow(result) ] <- "nobs"
    }
    coefs <- rrformat( result , digits=object@digits , width=object@width )
    print( coefs , quote=FALSE , justify="right" )
}
setMethod( "show" , "coeftab" , function(object) coeftab.show(object) )

coeftab <- function( ... , se=FALSE , se.inside=FALSE , nobs=TRUE , digits=2 , width=7 , rotate=FALSE , compare=FALSE ) {
    
    # se=TRUE outputs standard errors
    # se.inside=TRUE prints standard errors in parentheses in same column as estimates
    
    if ( se.inside==TRUE ) se <- TRUE
    
    # retrieve list of models
    L <- list(...)
    if ( is.list(L[[1]]) && length(L)==1 )
        L <- L[[1]]
    
    # retrieve model names from function call
    mnames <- match.call()
    mnames <- as.character(mnames)[2:(length(L)+1)]
    
    # count number of unique parameters
    param.names <- {}
    for ( i in 1:length(L) ) {
        c.names <- names( xcoef( L[[i]] ) )
        param.names <- unique( c( param.names , c.names ) )
    }
    # columns for standard errors
    if ( se==TRUE && se.inside==FALSE ) {
        for ( i in 1:length(L) ) {
            kse.names <- paste( names( xcoef( L[[i]] ) ) , ".se" , sep="" )
            param.names <- unique( c( param.names , kse.names ) )
        }
    }
    
    # make empty table
    nk <- length(param.names)
    d <- matrix( NA , ncol=nk )
    d <- data.frame(d)
    colnames(d) <- c( param.names )
    
    # loop over models and insert values
    for ( i in 1:length(L) ) {
        klist <- xcoef( L[[i]] )
        for ( j in 1:length(klist) ) {
            d[i,][ names( klist[j] ) ] <- as.numeric( round( klist[j] , digits ) )
        }
    }
    # insert standard errors
    if ( se==TRUE ) {
        for ( i in 1:length(L) ) {
            kse <- xse( L[[i]] )
            names(kse) <- names( xcoef( L[[i]] ) )
            for ( j in 1:length(kse) ) {
                if ( se.inside==FALSE )
                    # own column
                    d[i,][ paste(names(kse)[j],".se",sep="") ] <- as.numeric( round(kse[j],digits) )
                else
                    # combine with estimate
                    d[i,][ names(kse)[j] ] <- paste( formatC( as.real(d[i,][ names(kse)[j] ]) , digits=digits ) , " (" , formatC( as.real( kse[j] ) , digits=digits ) , ")" , sep="" )
            }
        }
    }
    
    # add model names to rows
    rownames(d) <- mnames
    
    # formatting for parenthetical standard errors
    if ( se.inside==TRUE && se==TRUE ) {
        comment(d) <- "Values in parentheses are quadratic estimate standard errors."
        colnames(d) <- paste( colnames(d) , "(se)" )
        for ( i in 1:nrow(d) ) {
            for ( j in 1:ncol(d) ) {
                d[i,j] <- ifelse( is.na(d[i,j]) , "" , d[i,j] )
            }
        }
    }
    
    # add nobs
    if ( nobs==TRUE ) {
        nobs <- sapply( L , xnobs )
    } else {
        nobs <- 0
    }
    
    # add AICc and weights
    if ( compare==TRUE ) {
        require(bbmle)
        d$AICc <- sapply( 1:length(L) , function(z) AICc( L[[z]] , nobs=d$nobs[z] ) )
        deltAICc <- d$AICc - min(d$AICc) 
        wAIC <- exp( -0.5*deltAICc )
        d$weight <- wAIC / sum( wAIC )
    }
    
    # return table
    if ( rotate==FALSE ) d <- t(d) # models along top is default
    new( "coeftab" , coefs=d , nobs=nobs , digits=digits , width=width )
}
