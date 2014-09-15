# sim -- simulates observations from fit map and map2stan models

setGeneric("sim",
function( fit , data , n=1000 , ... ) {
    predict(fit)
}
)

setMethod("sim", "map",
function( fit , data , n=1000 , post , ll=FALSE , refresh=0.1 , ... ) {
    # when ll=FALSE, simulates sampling, one sample for each sample in posterior
    # when ll=TRUE, computes loglik of each observation, for each sample in posterior
    
    ########################################
    # check arguments
    if ( class(fit)!="map" ) stop("Requires map fit")
    if ( missing(data) ) {
        data <- fit@data
    } else {
        # make sure all variables same length
        # weird vectorization errors otherwise
        #data <- as.data.frame(data)
    }
    
    if ( missing(post) ) 
        post <- extract.samples(fit,n=n)
    else
        n <- length(post[[1]])
    
    # get linear model values from link
    # use our posterior samples, so later parameters have right correlation structure
    # don't flatten result, so we end up with a named list, even if only one element
    pred <- link( fit , data=data , n=n , post=post , flatten=FALSE , refresh=refresh , ... )
    
    # extract likelihood, assuming it is first element of formula
    lik <- flist_untag(fit@formula)[[1]]
    # discover outcome
    outcome <- as.character(lik[[2]])
    # discover likelihood function
    flik <- as.character(lik[[3]][[1]])
    # get simulation partner function
    rlik <- flik
    if ( ll==FALSE ) substr( rlik , 1 , 1 ) <- "r"
    
    # check for aggreagted binomial, but only when ll==TRUE
    aggregated_binomial <- FALSE
    if ( flik=="dbinom" & ll==TRUE ) {
        ftemp <- flist_untag(fit@formula)[[1]]
        if ( class(ftemp[[3]][[2]])=="name" ) {
            aggregated_binomial <- TRUE
        }
    }
    
    # pull out parameters in likelihood
    pars <- vector(mode="list",length=length(lik[[3]])-1)
    for ( i in 1:length(pars) ) {
        pars[[i]] <- lik[[3]][[i+1]]
    }
    if ( aggregated_binomial==TRUE ) {
        # swap 'size' for '1' and outcome for vector of ones
        # will expand to correct number of cases later
        pars[[1]] <- 1
        data[['ones__']] <- rep( 1 , length(data[[outcome]]) )
        if ( refresh > 0 )
            message("Aggregated binomial counts detected. Splitting to 0/1 outcome for WAIC calculation.")
    }
    pars <- paste( pars , collapse=" , " )
    # build expression to evaluate
    n_cases <- length(data[[outcome]])
    if ( ll==FALSE ) {
        xeval <- paste( rlik , "(" , n_cases , "," , pars , ")" , collapse="" )
    } else {
        use_outcome <- outcome
        if ( aggregated_binomial==TRUE ) use_outcome <- 'ones__'
        xeval <- paste( rlik , "(" , use_outcome , "," , pars , ",log=TRUE )" , collapse="" )
    }
    
    # simulate outcomes
    sim_out <- matrix(NA,nrow=n,ncol=n_cases)
    for ( s in 1:n ) {
        # build environment
        elm <- pred[[1]][s,] # extract s-th sample case calculations
        elm <- list(elm)
        names(elm)[1] <- names(pred)[1]
        init <- list() # holds one row of samples across all params
        for ( j in 1:length(post) ) {
            par_name <- names(post)[ j ]
            dims <- dim( post[[par_name]] )
            # scalar
            if ( length(dims)<2 ) init[[par_name]] <- post[[par_name]][s]
            # vector
            if ( length(dims)==2 ) init[[par_name]] <- post[[par_name]][s,]
            # matrix
            if ( length(dims)==3 ) init[[par_name]] <- post[[par_name]][s,,]
        }#j
        e <- list( as.list(data) , as.list(init) , as.list(elm) )
        e <- unlist( e , recursive=FALSE )
        # evaluate
        sim_out[s,] <- eval(parse(text=xeval),envir=e)
    }
    
    # check for aggregated binomial outcome
    if ( aggregated_binomial==TRUE ) {
        # aggregated binomial with data for 'size'
        # need to split binomial counts in outcome into series of 0/1 outcomes
        # (1) sum 'size' variable in order to get number of new cases
        size_var <- data[[ as.character(ftemp[[3]][[2]]) ]]
        n_newcases <- sum(size_var)
        # (2) make new sim_out with expanded dimension
        sim_out_new <- matrix(NA,nrow=n,ncol=n_newcases)
        # (3) loop through each aggregated case and fill sim_out_new with bernoulli loglik
        current_newcase <- 1
        outcome_var <- data[[ outcome ]]
        for ( i in 1:n_cases ) {
            num_ones <- outcome_var[i]
            num_zeros <- size_var[i] - num_ones
            ll1 <- sim_out[,i]
            ll0 <- log( 1 - exp(ll1) )
            for ( j in 1:num_ones ) {
                sim_out_new[,current_newcase] <- ll1
                current_newcase <- current_newcase + 1
            }
            for ( j in 1:num_zeros ) {
                sim_out_new[,current_newcase] <- ll0
                current_newcase <- current_newcase + 1
            }
        }#i
        sim_out <- sim_out_new
    }
    
    # result
    return(sim_out)
}
)

setMethod("sim", "map2stan",
function( fit , data , n=0 , post , ... ) {
    
    ########################################
    # check arguments
    if ( missing(data) ) {
        data <- fit@data
    } else {
        # make sure all variables same length
        # weird vectorization errors otherwise
        #data <- as.data.frame(data)
    }
    
    if ( missing(post) ) 
        post <- extract.samples(fit,n=n)
    
    n <- dim(post[[1]])[1] # in case n=0 was passed
    
    # get linear model values from link
    # use our posterior samples, so later parameters have right correlation structure
    # don't flatten result, so we end up with a named list, even if only one element
    pred <- link( fit , data=data , n=n , post=post , flatten=FALSE , ... )
    
    # extract likelihood, assuming it is first element of formula
    lik <- flist_untag(fit@formula)[[1]]
    # discover outcome
    outcome <- as.character(lik[[2]])
    # discover likelihood function
    flik <- as.character(lik[[3]][[1]])
    # get simulation partner function
    rlik <- flik
    substr( rlik , 1 , 1 ) <- "r"
    # pull out parameters in likelihood
    pars <- vector(mode="list",length=length(lik[[3]])-1)
    for ( i in 1:length(pars) ) {
        pars[[i]] <- lik[[3]][[i+1]]
    }
    pars <- paste( pars , collapse=" , " )
    # build expression to evaluate
    n_cases <- length(data[[1]])
    xeval <- paste( rlik , "(" , n_cases , "," , pars , ")" , collapse="" )
    
    # simulate outcomes
    sim_out <- matrix( NA , nrow=n , ncol=n_cases )
    # need to select out s-th sample for each parameter in post
    # only need top-level parameters, so can coerce to data.frame?
    post <- as.data.frame(post)
    for ( s in 1:n ) {
        # build environment
        elm <- pred[[1]][s,] # extract s-th sample case calculations
        elm <- list(elm)
        names(elm)[1] <- names(pred)[1]
        e <- list( as.list(data) , as.list(post[s,]) , as.list(elm) )
        e <- unlist( e , recursive=FALSE )
        # evaluate
        sim_out[s,] <- eval(parse(text=xeval),envir=e)
    }
    
    return(sim_out)
}
)

postcheck <- function( fit , prob=0.9 , window=20 , ... ) {
    
    undot <- function( astring ) {
        astring <- gsub( "." , "_" , astring , fixed=TRUE )
        astring
    }
    
    pred <- link(fit)
    sims <- sim(fit)
    
    if ( class(pred)=="list" )
        if ( length(pred)>1 ) pred <- pred[[1]]
    
    # get outcome variable
    lik <- flist_untag(fit@formula)[[1]]
    outcome <- as.character(lik[[2]])
    if ( class(fit)=="map2stan" ) outcome <- undot(outcome)
    y <- fit@data[[ outcome ]]
    
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
    
    start <- 1
    end <- cases_per_page
    
    for ( i in 1:num_pages ) {
        
        end <- min( ny , end )
        window <- start:end
        
        set_nice_margins()
        plot( y[window] , xlab="case" , ylab=outcome , col=rangi2 , pch=16 , ylim=c( ymin , ymax ) , xaxt="n" )
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

# EXAMPLES/TEST
if ( FALSE ) {

library(rethinking)
data(cars)

fit <- map(
    alist(
        dist ~ dnorm(mu,sigma),
        mu <- a + b*speed
    ),
    data=cars,
    start=list(a=30,b=0,sigma=5)
)

pred <- link(fit,data=list(speed=0:30),flatten=TRUE)

sim.dist <- sim(fit,data=list(speed=0:30))

plot( dist ~ speed , cars )
mu <- apply( pred , 2 , mean )
mu.PI <- apply( pred , 2, PI )
y.PI <- apply( sim.dist , 2 , PI )
lines( 0:30 , mu )
shade( mu.PI , 0:30 )
shade( y.PI , 0:30 )

# now map2stan

cars$y <- cars$dist # Stan doesn't allow variables named 'dist'

m <- map2stan(
    alist(
        y ~ dnorm(mu,sigma),
        mu <- a + b*speed,
        a ~ dnorm(0,100),
        b ~ dnorm(0,10),
        sigma ~ dunif(0,100)
    ),
    data=cars,
    start=list(
        a=mean(cars$y),
        b=0,
        sigma=1
    ),
    warmup=500,iter=2500
)

pred <- link(m,data=list(speed=0:30),flatten=TRUE)

sim.dist <- sim(m,data=list(speed=0:30),n=0)

plot( dist ~ speed , cars )
mu <- apply( pred , 2 , mean )
mu.PI <- apply( pred , 2, PI )
y.PI <- apply( sim.dist , 2 , PI )
lines( 0:30 , mu )
shade( mu.PI , 0:30 )
shade( y.PI , 0:30 )


} #EXAMPLES
