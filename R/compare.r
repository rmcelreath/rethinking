# compare

# compare class definition and show method
setClass( "compareIC" , representation( output="data.frame" , dSE="matrix" ) )

compare.show <- function( object ) {
    r <- format_show( object@output , digits=c('default__'=1,'weight'=2,'SE'=2,'dSE'=2) )
    print( r )
}
setMethod( "show" , "compareIC" , function(object) compare.show(object) )

# new compare function, defaulting to WAIC
compare <- function( ... , n=1e3 , sort="WAIC" , WAIC=TRUE , refresh=0 ) {
    # retrieve list of models
    L <- list(...)
    if ( is.list(L[[1]]) && length(L)==1 )
        L <- L[[1]]
    
    # retrieve model names from function call
    mnames <- match.call()
    mnames <- as.character(mnames)[2:(length(L)+1)]
    
    dSE.matrix <- matrix( NA , nrow=length(L) , ncol=length(L) )
    if ( WAIC==FALSE ) {
        DIC.list <- lapply( L , function(z) DIC( z , n=n ) )
        pD.list <- sapply( DIC.list , function(x) attr(x,"pD") )
    } else {
        # use WAIC instead of DIC
        # WAIC is processed pointwise, so can compute SE of differences, and summed later
        WAIC.list <- lapply( L , function(z) WAIC( z , n=n , refresh=refresh , pointwise=TRUE ) )
        pD.list <- sapply( WAIC.list , function(x) sum(attr(x,"pWAIC")) )
        se.list <- sapply( WAIC.list , function(x) attr(x,"se") )
        DIC.list <- sapply( WAIC.list , sum )
        # compute SE of differences between adjacent models from top to bottom in ranking
        colnames(dSE.matrix) <- mnames
        rownames(dSE.matrix) <- mnames
        for ( i in 1:(length(L)-1) ) {
            for ( j in (i+1):length(L) ) {
                waic_ptw1 <- WAIC.list[[i]]
                waic_ptw2 <- WAIC.list[[j]]
                dSE.matrix[i,j] <- as.numeric( sqrt( length(waic_ptw1)*var( waic_ptw1 - waic_ptw2 ) ) )
                dSE.matrix[j,i] <- dSE.matrix[i,j]
            }#j
        }#i
    }
    
    DIC.list <- unlist(DIC.list)
    
    dDIC <- DIC.list - min( DIC.list )
    w.DIC <- ICweights( DIC.list )
    
    if ( WAIC==FALSE )
        result <- data.frame( DIC=DIC.list , pD=pD.list , dDIC=dDIC , weight=w.DIC )
    else {
        # find out which model has dWAIC==0
        topm <- which( dDIC==0 )
        dSEcol <- dSE.matrix[,topm]
        result <- data.frame( WAIC=DIC.list , pWAIC=pD.list , dWAIC=dDIC , 
                              weight=w.DIC , SE=se.list , dSE=dSEcol )
    }
    
    rownames(result) <- mnames
    
    if ( !is.null(sort) ) {
        if ( sort!=FALSE ) {
            if ( WAIC==FALSE & sort=="WAIC" ) sort <- "DIC"
            result <- result[ order( result[[sort]] ) , ]
        }
    }
    
    new( "compareIC" , output=result , dSE=dSE.matrix )
}

# plot method for compareIC results shows deviance in and expected deviance out of sample, for each model, ordered top-to-bottom by rank
setMethod("plot" , "compareIC" , function(x,y,xlim,SE=TRUE,dSE=TRUE,weights=FALSE,...) {
    dev_in <- x@output[[1]] - x@output[[2]]*2
    dev_out <- x@output[[1]]
    if ( !is.null(x@output[['SE']]) ) devSE <- x@output[['SE']]
    dev_out_lower <- dev_out - devSE
    dev_out_upper <- dev_out + devSE
    if ( weights==TRUE ) {
        dev_in <- ICweights(dev_in)
        dev_out <- ICweights(dev_out)
        dev_out_lower <- ICweights(dev_out_lower)
        dev_out_upper <- ICweights(dev_out_upper)
    }
    n <- length(dev_in)
    if ( missing(xlim) ) {
        xlim <- c(min(dev_in),max(dev_out))
        if ( SE==TRUE & !is.null(x@output[['SE']]) ) {
            xlim <- c(min(dev_in),max(dev_out_upper))
        }
    }
    main <- colnames(x@output)[1]
    set_nice_margins()
    dotchart( dev_in[n:1] , labels=rownames(x@output)[n:1] , xlab="deviance" , pch=16 , xlim=xlim , ... )
    points( dev_out[n:1] , 1:n )
    mtext(main)
    # standard errors
    if ( !is.null(x@output[['SE']]) & SE==TRUE ) {
        for ( i in 1:n ) {
            lines( c(dev_out_lower[i],dev_out_upper[i]) , rep(n+1-i,2) , lwd=0.75 )
        }
    }
    if ( !all(is.na(x@dSE)) & dSE==TRUE ) {
        # plot differences and stderr of differences
        dcol <- col.alpha("black",0.5)
        abline( v=dev_out[1] , lwd=0.5 , col=dcol )
        diff_dev_lower <- dev_out - x@output$dSE
        diff_dev_upper <- dev_out + x@output$dSE
        if ( weights==TRUE ) {
            diff_dev_lower <- ICweights(diff_dev_lower)
            diff_dev_upper <- ICweights(diff_dev_upper)
        }
        for ( i in 2:n ) {
            points( dev_out[i] , n+2-i-0.5 , cex=0.5 , pch=2 , col=dcol )
            lines( c(diff_dev_lower[i],diff_dev_upper[i]) , rep(n+2-i-0.5,2) , lwd=0.5 , col=dcol )
        }
    }
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
