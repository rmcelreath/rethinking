# compare

# compare class definition and show method
setClass( "compareIC" , slots=c( dSE="matrix" ) , contains="data.frame" )

#compare.show <- function( object ) {
#    r <- format_show( object@output , #digits=c('default__'=1,'weight'=2,'SE'=2,'dSE'=2) )
#    print( r )
#}
setMethod( "show" , "compareIC" , function(object) {
    r <- format_show( object , 
                      digits=c('default__'=1,'weight'=2,'SE'=2,'dSE'=2) )
    print( r )
} )

# new compare function, defaulting to WAIC
compare <- function( ... , n=1e3 , sort="WAIC" , func=WAIC , WAIC=TRUE , refresh=0 , warn=TRUE , result_order=c(1,5,3,6,2,4) , log_lik="log_lik" ) {

    # retrieve list of models
    L <- list(...)
    if ( is.list(L[[1]]) && length(L)==1 )
        L <- L[[1]]

    if ( length(L)==1 ) stop( "Need more than one model to compare." )
    
    # retrieve model names from function call
    mnames <- match.call()
    mnames <- as.character(mnames)[2:(length(L)+1)]

    # use substitute to deparse the func argument
    the_func <- func
    if ( class(the_func) != "character" )
        the_func <- deparse(substitute(func))
    
    # check class of fit models and warn when more than one class represented
    classes <- as.character(sapply( L , class ))
    if ( any(classes!=classes[1]) & warn==TRUE ) {
        warning("Not all model fits of same class.\nThis is usually a bad idea, because it implies they were fit by different algorithms.\nCheck yourself, before you wreck yourself.")
    }
    
    # check nobs for all models
    # if different, warn
    nobs_list <- try( sapply( L , nobs ) )
    if ( any(nobs_list != nobs_list[1]) & warn==TRUE ) {
        nobs_out <- paste( mnames , nobs_list , "\n" )
        nobs_out <- concat(nobs_out)
        warning(concat(
            "Different numbers of observations found for at least two models.\nModel comparison is valid only for models fit to exactly the same observations.\nNumber of observations for each model:\n",nobs_out))
    }
    
    dSE.matrix <- matrix( NA , nrow=length(L) , ncol=length(L) )
    # deprecate WAIC==TRUE/FALSE flag
    # catch it and convert to func intent
    if ( WAIC==FALSE ) func <- DIC # assume old code that wants DIC

    if ( the_func=="DIC" ) {
        IC.list <- lapply( L , function(z) DIC( z , n=n ) )
        p.list <- sapply( IC.list , function(x) attr(x,"pD") )
    }

    if ( the_func %in% c("WAIC","PSIS","LOO") ) {
        IC.list.pw <- lapply( L , function(z) do.call( the_func , list( z , n=n , refresh=refresh , pointwise=TRUE , warn=warn , log_lik=log_lik ) ) )
        p.list <- sapply( IC.list.pw , function(x) sum( x$penalty ) )
        se.list <- sapply( IC.list.pw , function(x) x$std_err[1] )
        IC.list <- sapply( IC.list.pw , function(x) sum( x[[1]] ) )
        # compute SE of differences between adjacent models from top to bottom in ranking
        colnames(dSE.matrix) <- mnames
        rownames(dSE.matrix) <- mnames
        for ( i in 1:(length(L)-1) ) {
            for ( j in (i+1):length(L) ) {
                ic_ptw1 <- IC.list.pw[[i]][[1]]
                ic_ptw2 <- IC.list.pw[[j]][[1]]
                dSE.matrix[i,j] <- as.numeric( sqrt( length(ic_ptw1)*var( ic_ptw1 - ic_ptw2 ) ) )
                dSE.matrix[j,i] <- dSE.matrix[i,j]
            }#j
        }#i
    }

    if ( !(the_func %in% c("DIC","WAIC","LOO","PSIS")) ) {
        # unrecognized IC function; just wing it
        IC.list <- lapply( L , function(z) func( z ) )
    }
    
    IC.list <- unlist(IC.list)
    
    dIC <- IC.list - min( IC.list )
    w.IC <- ICweights( IC.list )
    
    if ( the_func=="DIC" )
        result <- data.frame( DIC=IC.list , pD=p.list , dDIC=dIC , weight=w.IC )
    
    if ( the_func=="WAIC" ) {
        # find out which model has dWAIC==0
        topm <- which( dIC==0 )
        dSEcol <- dSE.matrix[,topm]
        result <- data.frame( WAIC=IC.list , pWAIC=p.list , dWAIC=dIC , 
                              weight=w.IC , SE=se.list , dSE=dSEcol )
    }
    if ( the_func=="LOO" | the_func=="PSIS" ) {
        topm <- which( dIC==0 )
        dSEcol <- dSE.matrix[,topm]
        result <- data.frame( PSIS=IC.list , pPSIS=p.list , dPSIS=dIC , 
                              weight=w.IC , SE=se.list , dSE=dSEcol )
    }
    
    if ( !(the_func %in% c("DIC","WAIC","LOO","PSIS")) ) {
        result <- data.frame( IC=IC.list , dIC=dIC , weight=w.IC )
    }
    
    rownames(result) <- mnames
    
    if ( !is.null(sort) ) {
        if ( sort!=FALSE ) {
            if ( the_func=="LOO" ) the_func <- "PSIS" # must match header
            if ( sort=="WAIC" ) sort <- the_func
            result <- result[ order( result[[sort]] ) , ]
        }
    }

    result <- result[ , result_order ]
    
    new( "compareIC" , result , dSE=dSE.matrix )
}

# plot method for compareIC results shows deviance in and expected deviance out of sample, for each model, ordered top-to-bottom by rank
setMethod("plot" , "compareIC" , function(x,y,xlim,SE=TRUE,dSE=TRUE,weights=FALSE,...) {
    dev_in <- x[[1]] - x[[5]]*2 # criterion - penalty*2
    dev_out <- x[[1]]
    if ( !is.null(x[['SE']]) ) devSE <- x[['SE']]
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
        if ( SE==TRUE & !is.null(x[['SE']]) ) {
            xlim[1] <- min(dev_in,dev_out_lower)
            xlim[2] <- max(dev_out_upper)
        }
    }
    main <- colnames(x)[1]
    set_nice_margins()
    dotchart( dev_in[n:1] , labels=rownames(x)[n:1] , xlab="deviance" , pch=16 , xlim=xlim , ... )
    points( dev_out[n:1] , 1:n )
    mtext(main)
    # standard errors
    if ( !is.null(x[['SE']]) & SE==TRUE ) {
        for ( i in 1:n ) {
            lines( c(dev_out_lower[i],dev_out_upper[i]) , rep(n+1-i,2) , lwd=0.75 )
        }
    }
    if ( !all(is.na(x@dSE)) & dSE==TRUE ) {
        # plot differences and stderr of differences
        dcol <- col.alpha("black",0.5)
        abline( v=dev_out[1] , lwd=0.5 , col=dcol )
        diff_dev_lower <- dev_out - x$dSE
        diff_dev_upper <- dev_out + x$dSE
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

if ( FALSE ) {
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
}#FALSE

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
    pulled_left = chimpanzees$pulled_left ,
    prosoc_left = chimpanzees$prosoc_left ,
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

# now ulam

m0 <- ulam(
    alist(
        pulled_left ~ dbinom(1,theta),
        logit(theta) <- a,
        a ~ dnorm(0,1)
    ) ,
    data=d , log_lik=TRUE )

m1 <- ulam(
    alist(
        pulled_left ~ dbinom(1,theta),
        logit(theta) <- a + bp*prosoc_left,
        a ~ dnorm(0,1),
        bp ~ dnorm(0,1)
    ) ,
    data=d , log_lik=TRUE )

m2 <- ulam(
    alist(
        pulled_left ~ dbinom(1,theta),
        logit(theta) <- a + bp*prosoc_left + bpc*condition*prosoc_left,
        a ~ dnorm(0,1),
        bp ~ dnorm(0,1),
        bpc ~ dnorm(0,1)
    ) ,
    data=d , log_lik=TRUE )

m3 <- ulam(
    alist(
        pulled_left ~ dbinom(1,theta),
        logit(theta) <- a + bp*prosoc_left + bc*condition + bpc*condition*prosoc_left,
        a ~ dnorm(0,1),
        bp ~ dnorm(0,1),
        bc ~ dnorm(0,1),
        bpc ~ dnorm(0,1)
    ) ,
    data=d , log_lik=TRUE )

m4 <- ulam(
    alist(
        pulled_left ~ dbinom(1,theta),
        logit(theta) <- a + aj[actor] + bp*prosoc_left + bc*condition + bpc*condition*prosoc_left,
        a ~ dnorm(0,1),
        aj[actor] ~ dnorm(0,sigma_actor),
        bp ~ dnorm(0,1),
        bc ~ dnorm(0,1),
        bpc ~ dnorm(0,1),
        sigma_actor ~ dcauchy(0,1)
    ) ,
    data=d , log_lik=TRUE )

( x1 <- compare(m0,m1,m2,m3,m4,log_lik="log_lik") )

plot(x1)


}
