# extracts n_divergent from stan fit
divergent <- function( fit , warmup=FALSE ) {
    x <- attr(fit,"cstanfit")$sampler_diagnostics(inc_warmup=warmup)
    sum(x[,,2])
}

# all diagnostics from Stan
# [1] "accept_stat__" "stepsize__"    "treedepth__"   "n_leapfrog__" 
# [5] "divergent__"   "energy__"  
#
# [1] "mean"    "se_mean" "sd"      "2.5%"    "25%"     "50%"     "75%"    
# [8] "97.5%"   "n_eff"   "Rhat" 
dashboard <- function( fit , warmup=FALSE , plot=TRUE , trank=TRUE ) {
    #if ( class(fit) %in% c("map2stan","ulam") ) fit <- fit@stanfit
    #x <- rstan::get_sampler_params(fit)
    x <- attr(fit,"cstanfit")$sampler_diagnostics(inc_warmup=warmup)
    n_chains <- dim(x)[2]
    y <- precis(fit,3) # get n_eff and Rhat
    n_samples <- dim(x)[1]*dim(x)[2]

    # display
    if ( plot==TRUE ) {
        set_nice_margins()
        par(cex.axis = 1 , cex.lab=1.2 )
        par(mfrow=c(2,2))

        # n_eff distribution
        Rhat_vals <- as.numeric( round( y[,'rhat'], 2 ) )
        plot( y[,'ess_bulk'] , (Rhat_vals) , xlab="number of effective samples" , ylab="Rhat" , ylim=c( 1 , max(1.1,(Rhat_vals),na.rm=TRUE) ) , log='y' )
        abline( v=0.1*n_samples , lty=1 , col="red" )
        abline( v=n_samples , lty=1 , col=grau() )
        abline( h=1 , lty=2 )

        # energy plot
        dens( x[,,3] , adj=0.1 , xlab="HMC energy" )
        mu <- mean(x[,,3])
        sig <- sd(x[,,3])
        curve( dnorm( x , mu , sig ) , add=TRUE , col=rangi2 , lwd=1.5 )

        # count of divergent iters
        n_divergent <- sum(x[,,2])
        plot( NULL , bty="n" , xlab=" " , ylab="" , xlim=c(0,1) , ylim=c(0,1) , xaxt="n" , yaxt="n" )
        text( 0.5 , 0.8 , n_divergent , cex=6 )
        text( 0.5 , 0.5 , "Divergent transitions" )
        if ( n_divergent==0 ) text( 0.5 , 0.2 , "Outlook good" , cex=1.5 )
        if ( n_divergent>0 & n_divergent<10 ) text( 0.5 , 0.2 , "Roll again" , cex=1.5 )
        if ( n_divergent>=10 ) text( 0.5 , 0.2 , "Check yourself before\nyou wreck yourself" , cex=1.5 )

        # three trace plots with lowest n_eff
        if ( trank==TRUE ) {
            trankplot( fit , pars="lp__" , lp=TRUE , add=TRUE )
        } else {
            traceplot( fit , pars="lp__" , lp=TRUE , add=TRUE )
        }
    }
    # invisible result
    invisible(x)
}
