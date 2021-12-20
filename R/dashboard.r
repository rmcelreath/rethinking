# extracts n_divergent from stan fit
divergent <- function( fit , warmup=FALSE ) {
    if ( class(fit) %in% c("map2stan","ulam") ) fit <- fit@stanfit
    x <- rstan::get_sampler_params(fit)
    if ( warmup==FALSE ) {
        nwarmup <- fit@stan_args[[1]]$warmup
        niter <- fit@stan_args[[1]]$iter
        n <- sapply( x , function(ch) sum(ch[(nwarmup+1):niter,"divergent__"]) )
    } else {
        n <- sapply( x , function(ch) sum(ch[,"divergent__"]) )
    }
    sum(n)
}

# all diagnostics from Stan
# [1] "accept_stat__" "stepsize__"    "treedepth__"   "n_leapfrog__" 
# [5] "divergent__"   "energy__"  
#
# [1] "mean"    "se_mean" "sd"      "2.5%"    "25%"     "50%"     "75%"    
# [8] "97.5%"   "n_eff"   "Rhat" 
dashboard <- function( fit , warmup=FALSE , plot=TRUE , trank=TRUE ) {
    if ( class(fit) %in% c("map2stan","ulam") ) fit <- fit@stanfit
    x <- rstan::get_sampler_params(fit)
    n_chains <- length(x)
    # remove warmup
    nwarmup <- fit@stan_args[[1]]$warmup
    niter <- fit@stan_args[[1]]$iter
    n_samples <- (niter - nwarmup)*n_chains
    if ( warmup==FALSE ) {
        for ( i in 1:n_chains ) {
            x[[i]] <- x[[i]][(nwarmup+1):niter,1:6]
        }
    }
    y <- summary(fit) # get n_eff and Rhat
    # merge chains
    if ( n_chains > 1 ) {
        x_temp <- x[[1]]
        for ( i in 2:n_chains )
            x_temp <- rbind(x_temp,x[[i]])
        x <- x_temp
    } else {
        x <- x[[1]]
    }

    # trace plot functions
    wstart <- floor( nwarmup * 0.5 )
    wend <- niter
    plot_make <- function(main, par, neff, ...) {
        ylim <- c(min(post[wstart:wend, , par]), max(post[wstart:wend, 
            , par]))
        plot(NULL, xlab = "", ylab = "", type = "l", xlim = c(wstart, 
            wend), ylim = ylim, ...)
        diff <- abs(ylim[1] - ylim[2])
        ylim <- ylim + c(-diff/2, diff/2)
        polygon(nwarmup * c(-1, 1, 1, -1), ylim[c(1, 1, 2, 2)], 
            col = grau(0.15), border = NA)
        neff_use <- neff
        mtext(paste("n_eff =", round(neff_use, 0)), 3, adj = 1, 
            cex = 0.9)
        mtext(main, 3, adj = 0, cex = 1)
    }
    plot_chain <- function(x, nc, ...) {
        lines(1:niter, x, col = col.alpha(rethink_palette[nc], 1), 
            lwd = 1)
    }

    pars <- "lp__" 
    # find pars with lowest n_eff
    #o <- order( y$summary[,9] )
    #pars <- c( pars , rownames(y$summary)[o[1:2]] )
    #print(pars)
    post <- extract(fit, pars = pars, permuted = FALSE, inc_warmup = TRUE)

    # display
    if ( plot==TRUE ) {
        set_nice_margins()
        par(cex.axis = 1 , cex.lab=1.2 )
        par(mfrow=c(2,2))

        # n_eff distribution
        Rhat_vals <- as.numeric( round( y$summary[,10], 2 ) )
        plot( y$summary[,9] , Rhat_vals , xlab="number of effective samples" , ylab="Rhat" , ylim=c( 0.995 , max(1.1,Rhat_vals,na.rm=TRUE) ) )
        abline( v=0.1*n_samples , lty=1 , col="red" )
        abline( v=n_samples , lty=1 , col=grau() )
        abline( h=1 , lty=2 )

        # energy plot
        dens( x[,6] , adj=0.1 , xlab="HMC energy" )
        mu <- mean(x[,6])
        sig <- sd(x[,6])
        curve( dnorm( x , mu , sig ) , add=TRUE , col=rangi2 , lwd=1.5 )

        # count of divergent iters
        plot( NULL , bty="n" , xlab=" " , ylab="" , xlim=c(0,1) , ylim=c(0,1) , xaxt="n" , yaxt="n" )
        text( 0.5 , 0.8 , sum(x[,4]) , cex=6 )
        text( 0.5 , 0.5 , "Divergent transitions" )
        if ( sum(x[,4])==0 ) text( 0.5 , 0.2 , "Outlook good" , cex=1.5 )
        if ( sum(x[,4])>10 ) text( 0.5 , 0.2 , "Check yourself before\nyou wreck yourself" , cex=1.5 )

        # three trace plots with lowest n_eff
        if ( trank==TRUE ) {
            trankplot( fit , pars="lp__" , lp=TRUE , add=TRUE )
        } else {
            plot_make( "log-probability" , "lp__" , y$summary["lp__","n_eff"] )
            for ( nc in 1:n_chains ) plot_chain( post[, nc , "lp__" ] , nc )
        }
    }
    # invisible result
    invisible(x)
}
