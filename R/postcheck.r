# graphical posterior validation checks for map and map2stan

postcheck <- function( fit , x , prob=0.89 , window=20 , n=1000 , col=rangi2 , ... ) {
    
    undot <- function( astring ) {
        astring <- gsub( "." , "_" , astring , fixed=TRUE )
        astring
    }
    
    pred <- link(fit,n=n)
    sims <- sim(fit,n=n)
    
    if ( class(pred)=="list" )
        if ( length(pred)>1 ) pred <- pred[[1]]
    
    # get outcome variable
    lik <- flist_untag(fit@formula)[[1]]
    dname <- as.character(lik[[3]][[1]])
    outcome <- as.character(lik[[2]])
    if ( class(fit)=="map2stan" ) outcome <- undot(outcome)
    y <- fit@data[[ outcome ]]
    
    # check for x-axis variable definition
    if ( !missing(x) ) {
        if ( x %in% names(fit@data) ) {
            x_var <- fit@data[[x]]
            r <- range(x_var)*c(0.9,1.1)
            x_seq <- seq(from=r[1],to=r[2],length.out=30)
        }
    }
    
    # compute posterior predictions for each case
    
    mu <- apply( pred , 2 , mean )
    mu.PI <- apply( pred , 2 , PI , prob=prob )
    y.PI <- apply( sims , 2 , PI , prob=prob )
    
    # figure out paging
    if ( length(window)==1 ) {
        window <- min(window,length(mu))
        num_pages <- ceiling( length(mu)/window )
        cases_per_page <- window
    }
    
    # display
    ny <- length(fit@data[[ outcome ]])
    ymin <- min(c(as.numeric(y.PI),mu,y))
    ymax <- max(c(as.numeric(y.PI),mu,y))
    
    # check for aggregated binomial context
    mumax <- max(c(as.numeric(mu.PI)))
    if ( ymax > 1 & mumax <= 1 & dname %in% c("dbinom","dbetabinom") ) {
        # probably aggregated binomial
        size_var <- as.character(lik[[3]][[2]])
        size_var <- fit@data[[ size_var ]]
        for ( i in 1:ny ) {
            y.PI[,i] <- y.PI[,i]/size_var[i]
            y[i] <- y[i]/size_var[i]
        }
        ymin <- 0
        ymax <- 1
    }
    if ( dname=="dordlogit" ) {
        # put mu and mu.PI on outcome scale
        ymin <- 1
        nlevels <- dim(pred)[3]
        mu <- sapply( 1:dim(pred)[2] , 
            function(i) { 
                temp <- t(pred[,i,])*1:7
                mean(apply(temp,2,sum))
            } )
        mu.PI <- sapply( 1:dim(pred)[2] , 
            function(i) { 
                temp <- t(pred[,i,])*1:7
                PI(apply(temp,2,sum),prob)
            } )
    }
    
    start <- 1
    end <- cases_per_page
    
    for ( i in 1:num_pages ) {
        
        end <- min( ny , end )
        window <- start:end
        
        set_nice_margins()
        plot( y[window] , xlab="case" , ylab=outcome , col=col , pch=16 , ylim=c( ymin , ymax ) , xaxt="n" , ... )
        axis( 1 , at=1:length(window) , labels=window )
        
        points( 1:length(window) , mu[window] )
        for ( x in 1:length(window) ) {
            lines( c(x,x) , mu.PI[,window[x]] )
            points( c(x,x) , y.PI[,window[x]] , cex=0.7 , pch=3 )
        }
        mtext( paste("Posterior validation check") )
        
        if ( num_pages > 1 ) {
            ask_old <- devAskNewPage(ask = TRUE)
            on.exit(devAskNewPage(ask = ask_old), add = TRUE)
        }
        
        start <- end + 1
        end <- end + cases_per_page
        
    }
    
    # invisible result
    result <- list(
        mean=mu,
        PI=mu.PI,
        outPI=y.PI
    )
    invisible( result )
}

