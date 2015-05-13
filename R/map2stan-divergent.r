# extracts n_divergent from stan fit
divergent <- function( fit , warmup=FALSE ) {
    if ( class(fit)=="map2stan" ) fit <- fit@stanfit
    x <- rstan::get_sampler_params(fit)
    if ( warmup==FALSE ) {
        nwarmup <- fit@stan_args[[1]]$warmup
        niter <- fit@stan_args[[1]]$iter
        n <- sapply( x , function(ch) sum(ch[(nwarmup+1):niter,5]) )
    } else {
        n <- sapply( x , function(ch) sum(ch[,5]) )
    }
    sum(n)
}

# all diagnostics from Stan
dashboard <- function( fit , warmup=FALSE ) {
    if ( class(fit)=="map2stan" ) fit <- fit@stanfit
    x <- rstan::get_sampler_params(fit)
    n_chains <- length(x)
    # remove warmup
    if ( warmup==FALSE ) {
        nwarmup <- fit@stan_args[[1]]$warmup
        niter <- fit@stan_args[[1]]$iter
        for ( i in 1:n_chains ) {
            x[[i]] <- x[[i]][(nwarmup+1):niter,1:5]
        }
    }
    # merge chains
    if ( n_chains > 1 ) {
        x_temp <- x[[1]]
        for ( i in 2:n_chains )
            x_temp <- rbind(x_temp,x[[i]])
        x <- x_temp
    } else {
        x <- x[[1]]
    }
    # display
    set_nice_margins()
    par(mfrow=c(2,2))
    # kernal of accept rate
    dens( x[,1] , xlab="accept rate" )
    # histogram of treedepth
    simplehist( x[,3] , xlab="tree depth" )
    # histogram of n_leapfrog
    simplehist( x[,4] , xlab="n leapfrog" )
    # count of divergent iters
    plot( NULL , bty="n" , xlab="n divergent" , ylab="" , xlim=c(0,1) , ylim=c(0,1) , xaxt="n" , yaxt="n" )
    text( 0.5 , 0.5 , sum(x[,5]) , cex=6 )
    # invisible result
    invisible(x)
}
