# method to implement ranked traceplots as described in Vehtari et al 2019 https://arxiv.org/abs/1903.08008
# given traces for a single parameter:
# rank draws across all chains
# then plot for each chain histogram of ranks

# convert matrix to a matrix of ranks (over entire matrix)
rank_mat <- function( x ) {
    if ( class(x)[1]=="numeric" ) x <- array( x , dim=c(length(x),1) )
    matrix( rank(x) , ncol=ncol(x) )
}

trankplot <- function( object , bins=30 , pars , chains , col=rethink_palette , alpha=1 , bg=col.alpha("black",0.15) , ask=TRUE , window , n_cols=3 , max_rows=5 , lwd=1.5 , lp=FALSE  , axes=FALSE , off=0 , add=FALSE , stacked=FALSE , ... ) {
    
    if ( !(class(object)[1] %in% c("map2stan","ulam","stanfit")) ) stop( "requires map2stan, ulam or stanfit object" )
    
    if ( class(object)[1] %in% c("map2stan","ulam") ) {
        object <- object@stanfit
    }

    # get all chains, not mixed, from stanfit
    # exclude warmup, because we'll rank only proper draws
    if ( missing(pars) )
        post <- extract(object,permuted=FALSE,inc_warmup=FALSE)
    else
        post <- extract(object,pars=pars,permuted=FALSE,inc_warmup=FALSE)
    
    # names
    dimnames <- attr(post,"dimnames")
    n_chains <- length(dimnames$chains)

    if ( n_chains==1 ) stop( "trankplot requires more than one chain." )

    if ( missing(chains) ) chains <- 1:n_chains
    n_chains <- length(chains)
    pars <- dimnames$parameters
    chain.cols <- rep_len(col,n_chains)
    # cut out "dev" and "lp__" and "log_lik"
    wdev <- which(pars=="dev")
    if ( length(wdev)>0 ) pars <- pars[-wdev]
    wlp <- which(pars=="lp__")
    if ( length(wlp)>0 & lp==FALSE ) pars <- pars[-wlp]
    wlp <- grep( "log_lik" , pars , fixed=TRUE )
    if ( length(wlp)>0 ) pars <- pars[-wlp]
    
    n_pars <- length( pars )

    # construct ranks
    # do this one parameter at a time
    ranks <- post
    n_samples <- dim(post)[1]
    for ( i in 1:n_pars ) {
        ranks[,,i] <- rank_mat( post[,, pars[i] ] )
    }
    breaks <- hist( ranks[,,1] , breaks=bins , plot=FALSE )$breaks
    h <- array( NA , dim=c( length(breaks)-1 , n_chains , n_pars ) )
    for ( i in 1:n_pars ) {
        for ( j in 1:n_chains ) {
            #print( hist( ranks[,j,i] , breaks=breaks , plot=FALSE )$counts )
            h[,j,i] <- hist( ranks[,j,i] , breaks=breaks , plot=FALSE )$counts
        }
    }

    # figure out grid and paging
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
    wstart <- 0
    wend <- max(breaks)
    if ( missing(window) ) window <- c(1,n_samples)
    if ( !missing(window) ) {
        wstart <- window[1]
        wend <- window[2]
    }
    
    # worker
    plot_make <- function( main , par , neff , ... ) {
        ylim <- range(h[,,par])
        if ( stacked==TRUE ) ylim[2] <- ylim[2] * ( length(chains) - 1 )
        if ( axes==TRUE )
            plot( NULL , xlab="" , ylab="" , bty="l" , xlim=range(breaks) , ylim=ylim , ... )
        else
            plot( NULL , xlab="" , ylab="" , bty="l" , xlim=range(breaks) , ylim=ylim , xaxt="n" , yaxt="n" , ... )
        neff_use <- neff[ names(neff)==main ]
        mtext( paste("n_eff =",round(neff_use,0)) , 3 , adj=1 , cex=0.9 )
        if ( main=="lp__" ) main <- "log-probability"
        mtext( main , 3 , adj=0 , cex=1 )
    }
    # make the trank
    nb <- length(breaks)
    plot_trank <- function( r , ... ) {
        # rank draws from all chains
        if ( stacked==FALSE ) {
            for ( i in chains ) {
                x <- c( breaks[1] , rep( breaks[2:(nb-1)] , each=2 ) , breaks[nb] )
                y <- rep( r[ 1:(nb-1) ,i] , each=2 )
                lines( x + (i-1)*off , y , col=col.alpha(chain.cols[i],alpha) , lwd=lwd )
            }#i
        } else {
            # stacked version
            ysum <- 0
            for ( i in chains ) {
                x <- c( breaks[1] , rep( breaks[2:(nb-1)] , each=2 ) , breaks[nb] )
                y <- rep( r[ 1:(nb-1) ,i] , each=2 )
                ysum <- y + ysum
                print(str(ysum))
                lines( x + (i-1)*off , ysum , col=col.alpha(chain.cols[i],alpha) , lwd=lwd )
            }#i
        }
    }
    
    # fetch n_eff
    n_eff <- summary(object)$summary[ , 'n_eff' ]
    
    # make window
    #set_nice_margins()
    if ( add==FALSE ) {
        par(mgp = c(0.5, 0.5, 0), mar = c(1.5, 1.5, 1.5, 1) + 0.1, tck = -0.02)
        par(mfrow=c(n_rows_per_page,n_cols))
    } 
    # draw traces
    n_ppp <- n_rows_per_page * n_cols # num pars per page
    for ( k in 1:n_pages ) {
        if ( k > 1 ) message( paste("Waiting to draw page",k,"of",n_pages) )
        for ( i in 1:n_ppp ) {
            pi <- i + (k-1)*n_ppp
            if ( pi <= n_pars ) {
                if ( pi == 2 ) {
                    if ( ask==TRUE & add==FALSE ) {
                        ask_old <- devAskNewPage(ask = TRUE)
                        on.exit(devAskNewPage(ask = ask_old), add = TRUE)
                    }
                }
                plot_make( pars[pi] , pi , n_eff , ... )
                plot_trank( h[ , , pi ] , ... )
                
            }
        }#i
        
    }#k
    
}
