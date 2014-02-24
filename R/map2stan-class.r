setClass("map2stan", representation( call = "language",
                                stanfit = "stanfit",
                                coef = "numeric",
                                vcov = "matrix",
                                data = "list",
                                start = "list",
                                pars = "character" ,
                                formula = "list" ,
                                formula_parsed = "list" ))

setMethod("coef", "map2stan", function(object) {
    object@coef
})

extract.samples <- function(object) {
    p <- extract(object@stanfit)
    # get rid of dev and lp__
    p[['dev']] <- NULL
    p[['lp__']] <- NULL
    # get rid of those ugly dimnames
    for ( i in 1:length(p) ) {
        attr(p[[i]],"dimnames") <- NULL
    }
    return(p)
}

# build ensemble of samples using DIC weights
ensemble <- function( ... ) {
    L <- list(...)
}

plotchains <- function(object , pars=names(object@start) , ...) {
    if ( class(object)=="map2stan" )
        rstan::traceplot( object@stanfit , ask=TRUE , pars=pars , ... )
}

plotpost <- function(object,n=1000,col=col.alpha("slateblue",0.3),cex=0.8,pch=16,...) {
    o <- as.data.frame(object)
    pairs(o[1:n,],col=col,cex=cex,pch=pch,...)
}

stancode <- function(object) {
    cat( object@stanfit@stanmodel@model_code )
    return( invisible( object@stanfit@stanmodel@model_code ) )
}

setMethod("vcov", "map2stan", function (object, ...) { object@vcov } )

setMethod("nobs", "map2stan", function (object, ...) { attr(object,"nobs") } )

setMethod("logLik", "map2stan",
function (object, ...)
{
    if(length(list(...)))
        warning("extra arguments discarded")
    val <- (-1)*attr(object,"deviance")/2
    attr(val, "df") <- length(object@coef)
    attr(val, "nobs") <- attr(object,"nobs")
    class(val) <- "logLik"
    val
  })
  
setMethod("deviance", "map2stan",
function (object, ...)
{
  attr(object,"deviance")
})

DIC <- function( object ) {
    val <- attr(object,"DIC")
    attr(val,"pD") <- attr(object,"pD")
    val
}

setMethod("show", "map2stan", function(object){

    cat("\nmap2stan model fit\n")
    
    cat("\nFormula:\n")
    for ( i in 1:length(object@formula) ) {
        print( object@formula[[i]] )
    }
    
    #cat("\nExpected values of fixed effects:\n")
    #print(coef(object))
    
    cat("\nLog-likelihood at expected values: ")
    cat(round(as.numeric(logLik(object)),2),"\n")
    
    cat("Deviance: ")
    cat(round(as.numeric(deviance(object)),2),"\n")
    
    cat("DIC: ")
    cat(round(as.numeric(DIC(object)),2),"\n")
    
    cat("Effective number of parameters (pD): ")
    cat(round(as.numeric(attr(object,"pD")),2),"\n")
    
  })

setMethod("summary", "map2stan", function(object){
    
    show(object@stanfit)
    
})

# resample from compiled map2stan fit
resample <- function( object , iter=1e4 , warmup=1000 , chains=1 , ... ) {
    if ( class(object)!="map2stan" )
        stop( "Requires map2stan fit or stanfit object" )
    
    init <- list()
    for ( i in 1:chains ) init[[i]] <- object@start
    fit <- stan( fit=object@stanfit , data=object@data , init=init , pars=object@pars , iter=iter , warmup=warmup , chains=chains , ... )
    result <- object
    result@stanfit <- fit
    return(result)
}

setMethod("plot" , "map2stan" , function(x,y,...) {
    #require(rstan)
    #rstan::traceplot( x@stanfit , ask=TRUE , pars=names(x@start) , ... )
    tracerplot(x,...)
})

setMethod("pairs" , "map2stan" , function(x, n=500 , alpha=0.7 , cex=0.7 , pch=16 , adj=1 , ...) {
    require(rstan)
    posterior <- extract.samples(x)
    panel.dens <- function(x, ...) {
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(usr[1:2], 0, 1.5) )
        h <- density(x,adj=adj)
        y <- h$y
        y <- y/max(y)
        abline( v=0 , col="gray" , lwd=0.5 )
        lines( h$x , y )
    }
    panel.2d <- function( x , y , ... ) {
        i <- sample( 1:length(x) , size=n )
        abline( v=0 , col="gray" , lwd=0.5 )
        abline( h=0 , col="gray" , lwd=0.5 )
        dcols <- densCols( x[i] , y[i] )
        dcols <- sapply( dcols , function(k) col.alpha(k,alpha) )
        points( x[i] , y[i] , col=dcols , ... )
    }
    panel.cor <- function( x , y , ... ) {
        k <- cor( x , y )
        cx <- sum(range(x))/2
        cy <- sum(range(y))/2
        text( cx , cy , round(k,2) , cex=2*exp(abs(k))/exp(1) )
    }
    pairs( posterior , cex=cex , pch=pch , upper.panel=panel.2d , lower.panel=panel.cor , diag.panel=panel.dens , ... )
})

# my trace plot function
tracerplot <- function( object , col=c("slateblue","black","orange","green") , alpha=0.7 , bg=gray(0.6,0.5) , ask=TRUE , ... ) {
    chain.cols <- col
    
    if ( class(object)!="map2stan" ) stop( "requires map2stan fit" )
    
    # get all chains, not mixed, from stanfit
    post <- extract(object@stanfit,permuted=FALSE,inc_warmup=TRUE)
    
    # names
    dimnames <- attr(post,"dimnames")
    chains <- dimnames$chains
    pars <- dimnames$parameters
    # cut out "dev" and "lp__"
    wdev <- which(pars=="dev")
    if ( length(wdev)>0 ) pars <- pars[-wdev]
    wlp <- which(pars=="lp__")
    if ( length(wdev)>0 ) pars <- pars[-wlp]
    
    # figure out grid and paging
    n_pars <- length( pars )
    n_cols=2
    n_rows=ceiling(n_pars/n_cols)
    n_rows_per_page <- n_rows
    paging <- FALSE
    n_pages <- 1
    if ( n_rows_per_page > 5 ) {
        n_rows_per_page <- 5
        n_pages <- ceiling(n_pars/(n_cols*n_rows_per_page))
        paging <- TRUE
    }
    n_iter <- object@stanfit@sim$iter
    n_warm <- object@stanfit@sim$warmup
    
    # worker
    plot_make <- function( main , par , neff , ... ) {
        ylim <- c( min(post[,,par]) , max(post[,,par]) )
        plot( NULL , xlab="sample" , ylab="position" , col=chain.cols[1] , type="l" , main=main , xlim=c(1,n_iter) , ylim=ylim , ... )
        # add polygon here for warmup region?
        diff <- abs(ylim[1]-ylim[2])
        ylim <- ylim + c( -diff/2 , diff/2 )
        polygon( n_warm*c(-1,1,1,-1) , ylim[c(1,1,2,2)] , col=bg , border=NA )
        mtext( paste("n_eff =",round(neff,0)) , 3 , adj=1 , cex=0.8 )
    }
    plot_chain <- function( x , nc , ... ) {
        lines( 1:n_iter , x , col=col.alpha(chain.cols[nc],alpha) )
    }
    
    # fetch n_eff
    n_eff <- summary(object@stanfit)$summary[,'n_eff']
    
    # make window
    mfrow_old <- par("mfrow")
    on.exit(par(mfrow = mfrow_old))
    par(mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1, 
            tck = -0.02)
    par(mfrow=c(n_rows_per_page,n_cols))
    
    # draw traces
    n_ppp <- n_rows_per_page * n_cols # num pars per page
    for ( k in 1:n_pages ) {
        if ( k > 1 ) message( paste("Waiting to draw page",k,"of",n_pages) )
        for ( i in 1:n_ppp ) {
            pi <- i + (k-1)*n_ppp
            if ( pi <= n_pars ) {
                if ( pi == 2 ) {
                    if ( ask==TRUE ) {
                        ask_old <- devAskNewPage(ask = TRUE)
                        on.exit(devAskNewPage(ask = ask_old), add = TRUE)
                    }
                }
                plot_make( pars[pi] , pi , n_eff[pi] , ... )
                for ( j in 1:length(chains) ) {
                    plot_chain( post[ , j , pi ] , j , ... )
                }#j
            }
        }#i
        
    }#k
    
}