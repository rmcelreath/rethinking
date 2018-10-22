setClass("map2stan", slots=c( call = "language",
                                model = "character",
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

stan_total_samples <- function(stanfit) {
    # find number of samples from stanfit
    iter <- stanfit@sim$iter
    warm <- stanfit@sim$warmup
    chains <- stanfit@sim$chains
    tot_samples <- (iter-warm)*chains
    return(tot_samples)
}

setMethod("extract.samples","map2stan",
function(object,n,...) {
    #require(rstan)
    p <- rstan::extract(object@stanfit,...)
    # get rid of dev and lp__
    p[['dev']] <- NULL
    p[['lp__']] <- NULL
    p[['log_lik']] <- NULL
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
    }
    return(p)
}
)

setMethod("extract.samples","stanfit",
function(object,...) {
    #require(rstan)
    p <- rstan::extract(object,...)
    # get rid of dev and lp__
    #p[['dev']] <- NULL
    #p[['lp__']] <- NULL
    # get rid of those ugly dimnames
    for ( i in 1:length(p) ) {
        attr(p[[i]],"dimnames") <- NULL
    }
    return(p)
}
)

plotchains <- function(object , pars=names(object@start) , ...) {
    if ( class(object)=="map2stan" )
        rstan::traceplot( object@stanfit , ask=TRUE , pars=pars , ... )
}

plotpost <- function(object,n=1000,col=col.alpha("slateblue",0.3),cex=0.8,pch=16,...) {
    o <- as.data.frame(object)
    pairs(o[1:n,],col=col,cex=cex,pch=pch,...)
}

setMethod("stancode", "map2stan",
function(object) {
    cat( object@model )
    return( invisible( object@model ) )
}
)
setMethod("stancode", "stanfit",
function(object) {
    cat( object@stanmodel@model_code )
    return( invisible( object@stanmodel@model_code ) )
}
)
setMethod("stancode", "list",
function(object) {
    cat( object$model )
    return( invisible( object$model ) )
}
)

setMethod("vcov", "map2stan", function (object, ...) { 
    #object@vcov 
    cov(as.data.frame(extract.samples(object,...)))
} )

setMethod("nobs", "map2stan", function (object, ...) { attr(object,"nobs") } )

setMethod("logLik", "map2stan",
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
  
setMethod("deviance", "map2stan",
function (object, ...)
{
    if ( is.null(attr(object,"deviance")) ) {
        return( as.numeric((-2)*logLik(object)) )
    } else {
        return( attr(object,"deviance") )
    }
})

stan_sampling_duration <- function(object) {
    if ( class(object)=="map2stan" | class(object)=="ulam" ) object <- object@stanfit
    dur <- get_elapsed_time(object)
    totals <- rep(NA,nrow(dur))
    for ( i in 1:nrow(dur) ) totals[i] <- sum( dur[i,] )
    dur <- cbind( dur , totals )
    colnames(dur)[3] <- "total"
    lab <- "seconds"
    if ( any(dur>60) ) {
        dur <- dur/60 # convert to minutes
        lab <- "minutes"
        if ( any(dur>60) ) {
            dur <- dur/60 # convert to hours
            lab <- "hours"
            if ( any(dur>24) ) {
                dur <- dur/24 # convert to days
                lab <- "days"
            }
        }
    }
    attr(dur,"units") <- lab
    return(dur)
}

setMethod("show", "map2stan", function(object){

    cat("map2stan model\n")
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
    
    #cat("\nExpected values of fixed effects:\n")
    #print(coef(object))
    
    if ( FALSE ) {
    cat("\nLog-likelihood at expected values: ")
    cat(round(as.numeric(logLik(object)),2),"\n")
    
    cat("Deviance: ")
    cat(round(as.numeric(deviance(object)),2),"\n")
    
    cat("DIC: ")
    cat(round(as.numeric(DIC(object)),2),"\n")
    
    cat("Effective number of parameters (pD): ")
    cat(round(as.numeric(attr(object,"pD")),2),"\n")
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

setMethod("summary", "map2stan", function(object){
    
    show(object@stanfit)
    
})

# resample from compiled map2stan fit
# can also run on multiple cores
resample <- function( object , ... ) {
    if ( class(object)!="map2stan" ) stop("Requires previous map2stan fit.")
    map2stan(object,...)
}
resample_old <- function( object , iter=1e4 , warmup=1000 , chains=1 , cores=1 , DIC=TRUE , WAIC=TRUE , rng_seed , data , ... ) {
    if ( !(class(object)%in%(c("map2stan"))) )
        stop( "Requires map2stan fit" )
    if ( missing(data) ) data <- object@data
    init <- list()
    if ( cores==1 | chains==1 ) {
        for ( i in 1:chains ) init[[i]] <- object@start
        fit <- stan( fit=object@stanfit , data=data , init=init , pars=object@pars , iter=iter , warmup=warmup , chains=chains , ... )
    } else {
        init[[1]] <- object@start
        #require(parallel)
        sys <- .Platform$OS.type
        if ( missing(rng_seed) ) rng_seed <- sample( 1:1e5 , 1 )
        if ( sys=='unix' ) {
            # Mac or Linux
            # hand off to mclapply
            sflist <- mclapply( 1:chains , mc.cores=cores ,
                function(chainid)
                    stan( fit=object@stanfit , data=data , init=init , pars=object@pars , iter=iter , warmup=warmup , chains=1 , seed=rng_seed, chain_id=chainid , ... )
            )
        } else {
            # Windows
            # so use parLapply instead
            CL = makeCluster(cores)
            fit <- object@stanfit
            #data <- object@data
            pars <- object@pars
            env0 <- list( fit=fit, data=data, pars=pars, rng_seed=rng_seed, iter=iter, warmup=warmup )
            clusterExport(cl = CL, c("iter","warmup","data", "fit", "pars", "rng_seed"), as.environment(env0))
            sflist <- parLapply(CL, 1:chains, fun = function(cid) {
                #require(rstan)
                stan(fit = fit, data = data, pars = pars, chains = 1, 
                  iter = iter, warmup = warmup, seed = rng_seed, 
                  chain_id = cid)
            })
        }
        # merge result
        fit <- sflist2stanfit(sflist)
    }
    
    result <- object
    result@stanfit <- fit
    
    # compute expected values of parameters
    s <- summary(fit)$summary
    s <- s[ -which( rownames(s)=="lp__" ) , ]
    s <- s[ -which( rownames(s)=="dev" ) , ]
    if ( !is.null(dim(s)) ) {
        coef <- s[,1]
        # compute variance-covariance matrix
        varcov <- matrix(NA,nrow=nrow(s),ncol=nrow(s))
        diag(varcov) <- s[,3]^2
    } else {
        coef <- s[1]
        varcov <- matrix( s[3]^2 , 1 , 1 )
        names(coef) <- names(result@start)
    }
    result@coef <- coef
    result@vcov <- varcov
    
    #DIC
    if ( DIC==TRUE ) {
        attr(result,"DIC") <- NULL
        dic_calc <- DIC(result)
        pD <- attr(dic_calc,"pD")
        attr(result,"DIC") <- dic_calc
        attr(result,"pD") <- pD
        attr(result,"deviance") <- dic_calc - 2*pD
    } else {
        # clear out any old DIC calculation
        attr(result,"DIC") <- NULL
        attr(result,"pD") <- NULL
        attr(result,"deviance") <- NULL
    }
    
    #WAIC
    if ( WAIC==TRUE ) {
        attr(result,"WAIC") <- NULL
        waic_calc <- try(WAIC(result,n=0))
        attr(result,"WAIC") <- waic_calc
    } else {
        # clear out any old WAIC calculation
        attr(result,"WAIC") <- NULL
    }
    
    return(result)
}

setMethod("plot" , "map2stan" , function(x,y,...) {
    #require(rstan)
    #rstan::traceplot( x@stanfit , ask=TRUE , pars=names(x@start) , ... )
    tracerplot(x,...)
})

setMethod("pairs" , "map2stan" , function(x, n=500 , alpha=0.7 , cex=0.7 , pch=16 , adj=1 , pars , ...) {
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
#rethink_palette <- c("#5BBCD6","#F98400","#F2AD00","#00A08A","#FF0000")
rethink_palette <- c("#8080FF","#F98400","#F2AD00","#00A08A","#FF0000")
rethink_cmyk <- c(col.alpha("black",0.25),"cyan")
tracerplot <- function( object , pars , col=rethink_palette , alpha=1 , bg=col.alpha("black",0.15) , ask=TRUE , window , n_cols=3 , max_rows=5 , lwd=0.5 , ... ) {
    
    if ( !(class(object) %in% c("map2stan","stanfit")) ) stop( "requires map2stan or stanfit fit object" )
    
    if ( class(object)=="map2stan" ) object <- object@stanfit

    # get all chains, not mixed, from stanfit
    if ( missing(pars) )
        post <- extract(object,permuted=FALSE,inc_warmup=TRUE)
    else
        post <- extract(object,pars=pars,permuted=FALSE,inc_warmup=TRUE)
    
    # names
    dimnames <- attr(post,"dimnames")
    chains <- dimnames$chains
    pars <- dimnames$parameters
    chain.cols <- rep_len(col,length(chains))
    # cut out "dev" and "lp__"
    wdev <- which(pars=="dev")
    if ( length(wdev)>0 ) pars <- pars[-wdev]
    wlp <- which(pars=="lp__")
    if ( length(wdev)>0 ) pars <- pars[-wlp]
    
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
    wstart <- 1
    wend <- n_iter
    if ( !missing(window) ) {
        wstart <- window[1]
        wend <- window[2]
    }
    
    # worker
    plot_make <- function( main , par , neff , ... ) {
        ylim <- c( min(post[wstart:wend,,par]) , max(post[wstart:wend,,par]) )
        plot( NULL , xlab="" , ylab="" , type="l" , xlim=c(wstart,wend) , ylim=ylim , ... )
        # add polygon here for warmup region?
        diff <- abs(ylim[1]-ylim[2])
        ylim <- ylim + c( -diff/2 , diff/2 )
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
                for ( j in 1:length(chains) ) {
                    plot_chain( post[ , j , pi ] , j , ... )
                }#j
            }
        }#i
        
    }#k
    
}

stanergy <- function( x , colscheme="blue" , binwidth=NULL , merge_chains=FALSE ) {
    library(bayesplot)
    if ( class(x)=="map2stan" ) x <- x@stanfit
    if ( class(x)!="stanfit" ) stop("needs a stanfit or map2stan object")
    np <- nuts_params(x)
    color_scheme_set(colscheme) # we hates the ggplot
    mcmc_nuts_energy( np , merge_chains = merge_chains , binwidth=binwidth )
}

# function to do mcmc_parcoord to help find divergences
divergence_tracker <- function( x , no_lp=TRUE , pars , ... ) {
    require(bayesplot)
    if ( class(x)=="map2stan" ) x <- x@stanfit
    np <- nuts_params(x)
    draws <- as.array(x)
    if ( missing(pars) ) {
        pars <- dimnames(x)$parameters
        if ( no_lp==TRUE ) {
            pars <- pars[ -which(pars=="lp__") ]
        }
    }
    bayesplot::mcmc_parcoord(draws,np=np,pars=pars,...)
}
# divergence_tracker( m5 )
