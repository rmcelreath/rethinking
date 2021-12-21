setClass("ulam", slots=c( call = "language",
                                model = "character",
                                #stanfit = "stanfit",
                                coef = "numeric",
                                vcov = "matrix",
                                data = "list",
                                start = "list",
                                pars = "character" ,
                                formula = "list" ,
                                formula_parsed = "list" ))

setMethod("coef", "ulam", function(object) {
    object@coef
})

setMethod("precis", "ulam",
function( object , depth=1 , pars , prob=0.89 , digits=2 , sort=NULL , decreasing=FALSE , lp__=FALSE , omit=NULL ,... ) {
    low <- (1-prob)/2
    upp <- 1-low

    # when fit with cmdstan, all parameters/variable in result
    # so want to filter at minimum by object@pars
    if ( missing(pars) ) pars <- object@pars

    if ( !is.null(attr(object,"stanfit")) ) {
        result <- summary( attr(object,"stanfit") ,pars=pars,probs=c(low,upp))$summary[,c(1,3:7)]
        result <- as.data.frame( result )
    }
    if ( !is.null(attr(object,"cstanfit")) ) {
        return( precis( attr(object,"cstanfit") , depth=depth, pars=pars , prob=prob, omit=omit , ... ) )
    }

    banlist <- c("dev","lp__")
    if ( lp__==TRUE ) banlist <- c("dev")
    idx <- which( rownames(result) %in% banlist )
    idx2 <- grep( "log_lik[" , rownames(result) , fixed=TRUE )
    if ( length(idx2)>0 ) idx <- c( idx , idx2 )
    if ( length(idx)>0 ) {
        # remove dev and lp__ and log_lik from table
        result <- result[ -idx , ]
    }
    # any pars to omit?
    if ( !is.null(omit) ) {
        for ( k in 1:length(omit) ) {
            idx <- grep( omit[k] , rownames(result) , fixed=TRUE )
            if ( length(idx)>0 ) result <- result[ -idx , ]
        }
    }

    result <- precis_format( result , depth , sort , decreasing )

    return( new( "precis" , result , digits=digits ) )
})

# models fit with cmdstan=TRUE include all parameters/variables
# so need to trim what is returned using object@pars
extract_post_ulam <- 
function(object,n,clean=TRUE,pars,...) {
    #require(rstan)
    if ( missing(pars) & clean==TRUE ) pars <- object@pars
    p <- rstan::extract(object@stanfit,pars=pars,...)
    # get rid of dev and lp__
    if ( clean==TRUE ) {
        p[['dev']] <- NULL
        p[['lp__']] <- NULL
        p[['log_lik']] <- NULL
    }
    # get rid of those ugly dimnames
    for ( i in 1:length(p) ) {
        attr(p[[i]],"dimnames") <- NULL
    }
    if ( !missing(n) ) {
        tot_samples <- stan_total_samples(object@stanfit)
        n <- min(n,tot_samples)
        for ( i in 1:length(p) ) {
            n_dims <- length( dim(p[[i]]) )
            if ( n_dims==1 ) p[[i]] <- p[[i]][1:n]
            if ( n_dims==2 ) p[[i]] <- p[[i]][1:n,]
            if ( n_dims==3 ) p[[i]] <- p[[i]][1:n,,]
        }
    } else {
        n <- stan_total_samples(object@stanfit)
    }

    model_name <- match.call()[[2]]
    attr(p,"source") <- concat( "ulam posterior: ", n , " samples from " , model_name )

    return(p)
}
setMethod("extract.samples","ulam",extract_post_ulam)

setMethod("stancode", "ulam",
function(object) {
    cat( object@model )
    return( invisible( object@model ) )
}
)

setMethod("vcov", "ulam", function (object, ...) { 
    #object@vcov 
    cov(as.data.frame(extract.samples(object,...)))
} )

setMethod("nobs", "ulam", function (object, ...) { attr(object,"nobs") } )

setMethod("logLik", "ulam",
function (object, ...)
{
    if(length(list(...)))
        warning("extra arguments discarded")
    if ( is.null(attr(object,"deviance") ) ) {
        val <- attr( WAIC(object) , "lppd" )
    } else {
        val <- (-1)*attr(object,"deviance")/2
    }
    attr(val, "df") <- length(object@coef)
    attr(val, "nobs") <- attr(object,"nobs")
    class(val) <- "logLik"
    val
  })
  
setMethod("deviance", "ulam",
function (object, ...)
{
    if ( is.null(attr(object,"deviance")) ) {
        return( as.numeric((-2)*logLik(object)) )
    } else {
        return( attr(object,"deviance") )
    }
})

setMethod("show", "ulam", function(object){

    cat("Hamiltonian Monte Carlo approximation\n")
    iter <- object@stanfit@sim$iter
    warm <- object@stanfit@sim$warmup
    chains <- object@stanfit@sim$chains
    chaintxt <- " chain\n"
    if ( chains>1 ) chaintxt <- " chains\n"
    tot_samples <- (iter-warm)*chains
    cat(concat( tot_samples , " samples from " , chains , chaintxt ))
    
    dur <- stan_sampling_duration(object)
    lab <- attr(dur,"units")
    attr(dur,"units") <- NULL
    cat(concat("\nSampling durations (",lab,"):\n"))
    print(round(dur,2))

    cat("\nFormula:\n")
    for ( i in 1:length(object@formula) ) {
        print( object@formula[[i]] )
    }
    
    if ( !is.null(attr(object,"WAIC")) ) {
        waic <- attr(object,"WAIC")
        use_waic <- sum(waic)
        cat("\nWAIC (SE): ")
        cat( concat(round(as.numeric(use_waic),2) , " (" , round(as.numeric(attr(waic,"se")),1) , ")" , "\n" ) )
        
        cat("pWAIC: ")
        use_pWAIC <- sum( unlist(attr(waic,"pWAIC")) )
        cat( round(as.numeric(use_pWAIC),2) , "\n" )
    }
    
  })

setMethod("summary", "ulam", function(object){
    
    show(object@stanfit)
    
})

setMethod("pairs" , "ulam" , function(x, n=200 , alpha=0.7 , cex=0.7 , pch=16 , adj=1 , pars , ...) {
    #require(rstan)
    if ( missing(pars) )
        posterior <- extract.samples(x)
    else
        posterior <- extract.samples(x,pars=pars)
    #if ( !missing(pars) ) {
    #    # select out named parameters
    #    p <- list()
    #    for ( k in pars ) p[[k]] <- posterior[[k]]
    #    posterior <- p
    #}
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
traceplot_ulam <- function( object , pars , chains , col=rethink_palette , alpha=1 , bg=col.alpha("black",0.15) , ask=TRUE , window , trim=100 , n_cols=3 , max_rows=5 , lwd=0.5 , lp=FALSE , ... ) {
    
    if ( !(class(object) %in% c("map2stan","ulam","stanfit")) ) stop( "requires map2stan or stanfit fit object" )
    
    if ( class(object) %in% c("map2stan","ulam") ) object <- object@stanfit

    # get all chains, not mixed, from stanfit
    if ( missing(pars) ) {
        post <- extract(object,permuted=FALSE,inc_warmup=TRUE)
        dimnames <- attr(post,"dimnames")
        pars <- dimnames$parameters
        # cut out "dev" and "lp__" and "log_lik"
        wdev <- which(pars=="dev")
        if ( length(wdev)>0 ) pars <- pars[-wdev]
        wlp <- which(pars=="lp__")
        if ( length(wlp)>0 & lp==FALSE ) pars <- pars[-wlp]
        wlp <- grep( "log_lik" , pars , fixed=TRUE )
        if ( length(wlp)>0 ) pars <- pars[-wlp]
    } else
        post <- extract(object,pars=pars,permuted=FALSE,inc_warmup=TRUE)
    
    # names
    dimnames <- attr(post,"dimnames")
    n_chains <- length(dimnames$chains)
    if ( missing(chains) ) chains <- 1:n_chains
    chain.cols <- rep_len(col,n_chains)
    
    # figure out grid and paging
    n_pars <- length( pars )
    n_rows=ceiling(n_pars/n_cols)
    n_rows_per_page <- n_rows
    paging <- FALSE
    n_pages <- 1
    if ( n_rows_per_page > max_rows ) {
        n_rows_per_page <- max_rows
        n_pages <- ceiling(n_pars/(n_cols*n_rows_per_page))
        paging <- TRUE
    }
    n_iter <- object@sim$iter
    n_warm <- object@sim$warmup
    n_samples_extracted <- dim( post )[1]
    wstart <- 1
    wend <- n_iter

    if ( !missing(window) ) {
        wstart <- window[1]
        wend <- window[2]
    }

    show_warmup <- TRUE
    if ( missing(window) ) {
        if ( n_iter > n_samples_extracted ) {
            # probably no warmup saved
            wend <- n_samples_extracted
            show_warmup <- FALSE
            trim <- 1 # no trim when warmup not shown
            n_iter <- n_samples_extracted
        }
        window <- c(trim,n_iter)
    }

    # worker
    plot_make <- function( main , par , neff , ... ) {
        ylim <- c( min(post[wstart:wend,,pars[par]]) , max(post[wstart:wend,,pars[par]]) )
        plot( NULL , xlab="" , ylab="" , type="l" , xlim=c(wstart,wend) , ylim=ylim , ... )
        # add polygon here for warmup region?
        diff <- abs(ylim[1]-ylim[2])
        ylim <- ylim + c( -diff/2 , diff/2 )
        if ( show_warmup==TRUE )
            polygon( n_warm*c(-1,1,1,-1) , ylim[c(1,1,2,2)] , col=bg , border=NA )
        neff_use <- neff[ names(neff)==main ]
        mtext( paste("n_eff =",round(neff_use,0)) , 3 , adj=1 , cex=0.9 )
        mtext( main , 3 , adj=0 , cex=1 )
    }
    plot_chain <- function( x , nc , ... ) {
        lines( 1:n_iter , x , col=col.alpha(chain.cols[nc],alpha) , lwd=lwd )
    }
    
    # fetch n_eff
    n_eff <- summary(object)$summary[ , 'n_eff' ]
    
    # make window
    #set_nice_margins()
    par(mgp = c(0.5, 0.5, 0), mar = c(1.5, 1.5, 1.5, 1) + 0.1, 
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
                plot_make( pars[pi] , pi , n_eff , ... )
                for ( j in 1:n_chains ) {
                    if ( j %in% chains )
                        plot_chain( post[ , j , pars[pi] ] , j , ... )
                }#j
            }
        }#i
        
    }#k
    
}
setMethod("traceplot", "ulam" , function(object,...) traceplot_ulam(object,...) )

setMethod( "plot" , "ulam" , function(x,depth=1,...) precis_plot(precis(x,depth=depth),...) )

setMethod("nobs", "ulam", function (object, ...) { 
    z <- attr(object,"nobs") 
    if ( is.null(z) ) {
        # try to get nobs from a link function
        link_funcs <- object@formula_parsed$link_funcs
        if ( length(link_funcs)>0 ) {
            z <- link_funcs[[1]]$N
        }
    }
    return(z)
} )



