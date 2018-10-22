# chapter 7

axis_unscale <- function( side=1 , at , orig , factor , ... ) {
    if ( missing(orig) )
        atx <- at / factor
    else
        atx <- ( at - mean(orig) ) / sd(orig)
    axis( side , at=atx , labels=at , ... )
}

brain_plot <- function( fit , atx=c(35,47,60) , aty=c(450,900,1300) , xlim , ylim , npts=100 ) {
    R2 <- R2_is_bad(fit)
    if ( is.nan(R2) ) R2 <- 1
    post <- extract.samples(fit)
    n <- dim( post$b )[2]
    if ( is.null(n) ) n <- 1

    if ( missing(xlim) ) 
        xlim <- range(d$mass_std)
    else
        xlim <- ( xlim - mean(d$mass) ) / sd(d$mass)
    if ( missing(ylim) ) 
        ylim <- range(d$brain_std)
    else
        ylim <- ylim / max(d$brain)

    plot( d$brain_std ~ d$mass_std , xaxt="n" , yaxt="n" , xlab="body mass (kg)" , ylab="brain volume (cc)" , col=rangi2 , pch=16 , xlim=xlim , ylim=ylim )
    axis_unscale( 1 , atx , d$mass )
    axis_unscale( 2 , at=aty , factor=max(d$brain) )
    xseq <- seq(from=xlim[1]-0.2,to=xlim[2]+0.2,length.out=npts)
    l <- link(fit,data=list(mass_std=xseq),refresh=0)
    mu <- apply(l,2,mean)
    lines( xseq , mu , lwd=2 )
    ci <- apply(l,2,PI)
    shade( ci , xseq )
    model_name <- deparse( match.call()[[2]] )
    mtext( concat(model_name,": R^2 = ",round(R2,2)) , adj=0 )
}

brain_loo_plot <- function( fit , atx=c(35,47,60) , aty=c(450,900,1300) , xlim , ylim , npts=100 ) {

    post <- extract.samples(fit)
    n <- dim( post$b )[2]
    if ( is.null(n) ) n <- 1

    if ( missing(xlim) ) 
        xlim <- range(d$mass_std)
    else
        xlim <- ( xlim - mean(d$mass) ) / sd(d$mass)
    if ( missing(ylim) ) 
        ylim <- range(d$brain_std)
    else
        ylim <- ylim / max(d$brain)

    plot( d$brain_std ~ d$mass_std , xaxt="n" , yaxt="n" , xlab="body mass (kg)" , ylab="brain volume (cc)" , col=rangi2 , pch=16 , xlim=xlim , ylim=ylim )
    axis_unscale( 1 , atx , d$mass )
    axis_unscale( 2 , at=aty , factor=max(d$brain) )

    d <- as.data.frame(fit@data)

    for ( i in 1:nrow(d) ) {
        di <- d[ -i , ]
        m_temp <- quap( fit@formula , data=di , start=list(b=rep(0,n)) )
        xseq <- seq(from=xlim[1]-0.2,to=xlim[2]+0.2,length.out=npts)
        l <- link(m_temp,data=list(mass_std=xseq),refresh=0)
        mu <- apply(l,2,mean)
        lines( xseq , mu , lwd=2 , col=col.alpha("black",0.3) )
    }

    model_name <- deparse( match.call()[[2]] )
    mtext( model_name , adj=0 )
}
