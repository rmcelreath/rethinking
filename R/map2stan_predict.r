# posterior predictions for map2stan model fits
# defaults to using original data

# to do:
# (*) fix dzipois vectorization hack --- need to fix dzipois itself
# (*) fix dordlogit in WAIC calculation -- not sure what issue is

# new link function that doesn't invoke Stan
setMethod("link", "map2stan",
function( fit , data , n=1000 , probs=NULL , refresh=0.1 , replace=list() , flatten=TRUE , ... ) {

    if ( class(fit)!="map2stan" ) stop("Requires map2stan fit")
    if ( missing(data) ) data <- fit@data
    
    post <- extract.samples(fit)
    
    # replace with any elements of replace list
    if ( length( replace ) > 0 ) {
        for ( i in 1:length(replace) ) {
            post[[ names(replace)[i] ]] <- replace[[i]]
        }
    }
    
    lm <- fit@formula_parsed$lm
    lik <- fit@formula_parsed$lik
    n_lm <- length(lm)
    n_lik <- length(lik)
    lm_names <- c()
    lm_lik <- c()
    if ( n_lm>0 ) {
        for ( i in 1:n_lm ) {
            lm_names <- c( lm_names , lm[[i]]$parameter )
            if ( n_lik > 0 ) {
                # find corresponding likelihood distribution
                for ( j in 1:n_lik ) {
                    if ( lik[[j]]$N_name == lm[[i]]$N_name ) {
                        lm_lik <- c( lm_lik , lik[[j]]$likelihood )
                    }
                }#j
            }
        }#i
    } else {
        stop( "There appear to be no linear models here" )
    }
    
    # number of samples
    n_samples <- dim( post[[1]] )[1]
    if ( n == 0 ) n <- n_samples # special flag for all samples in fit
    if ( n_samples < n ) n <- n_samples
    
    lm_out <- list()
    liks <- fit@formula_parsed$likelihood
    for ( i in 1:n_lm ) {
        # find out how many cases by finding corresponding likelihood and counting cases of outcome in data
        # can't just use old number of cases (from fit), because might have new user data
        n_cases <- 1
        K <- 1 # max level for ordered outcome
        N_name <- lm[[i]]$N_name
        outcome <- NA
        for ( j in 1:length(liks) ) {
            if ( liks[[j]]$N_name == N_name ) outcome <- liks[[j]]$outcome
        }
        if ( !is.na(outcome) ) {
            n_cases <- length( data[[ outcome ]] )
            if ( n_cases > 0 )
                K <- max( fit@data[[ outcome ]] ) # get max level from original data
            else {
                # outcome not found, so just fill with 1s
                n_cases <- length( data[[1]] ) # just hope is informative
                data[[outcome]] <- rep(1,n_cases)
            }
        }
        # empty array with right dims
        lm_out[[ lm_names[i] ]] <- array( 0 , dim=c( n , n_cases ) )
        
        if ( lm_lik[[i]]=="ordered_logistic" ) {
            # need another dimension in result, because outcome is a vector
            lm_out[[ lm_names[i] ]] <- array( 0 , dim=c( n , n_cases , K ) )
        }
    }
    
    init <- fit@start
    
    # loop over samples and compute each case for each linear model
    ref_inc <- floor(n*refresh)
    ref_next <- ref_inc
    
    # prep reused objects
    rhs <- list()
    for ( k in 1:n_lm ) {
        # ready linear model code
        rhs0 <- fit@formula_parsed$lm[[k]]$RHS
        rhs0 <- gsub( "[i]" , "" , rhs0 , fixed=TRUE )
        rhs[[k]] <- rhs0
    }
    
    for ( i in 1:n ) {
    
        # refresh progress display
        if ( refresh > 0 ) {
            if ( i == ref_next ) {
                msg <- paste( "[" , i , "/" , n , "]\r" , collapse=" " )
                cat(msg)
                ref_next <- i + ref_inc
                if ( ref_next > n ) ref_next <- n
            }
        }
        
        # build inits
        for ( j in 1:length(init) ) {
            par_name <- names(init)[ j ]
            dims <- dim( post[[par_name]] )
            # scalar
            if ( length(dims)==1 ) init[[par_name]] <- post[[par_name]][i]
            # vector
            if ( length(dims)==2 ) init[[par_name]] <- post[[par_name]][i,]
            # matrix
            if ( length(dims)==3 ) init[[par_name]] <- post[[par_name]][i,,]
        }#j
        
        # loop over linear models and compute by pushing samples through stanfit
        for ( k in 1:n_lm ) {
        
            # ready environment
            e <- list( as.list(data) , as.list(init) )
            e <- unlist( e , recursive=FALSE )
            
            # evaluate
            r <- eval(parse(text=rhs[[k]]),envir=e)
            flink <- fit@formula_parsed$lm[[k]]$link
            if ( flink=="log" ) r <- exp(r)
            if ( flink=="logit" ) r <- logistic(r)
            
            # special processing for ordered logit
            # this needs to be sped up
            if ( !is.null(lm_lik) ) {
                if ( lm_lik[k] == "ordered_logistic" ) {
                    # need vector of probabilities as result
                    cuts <- init[['cutpoints']]
                    v <- predict.ordlogit( r , cuts )
                    r <- t(v)
                }
            }
            
            # store
            if ( lm_lik[k] == "ordered_logistic" ) {
                lm_out[[ lm_names[k] ]][i,,] <- r
            } else {
                lm_out[[ lm_names[k] ]][i,] <- r
            }
        }#k
        
    }#i
    
    if ( !is.null(probs) ) {
        #NYI
    }
    
    if ( refresh>0 ) cat("\n")
    
    if ( flatten==TRUE )
        if ( length(lm_out)==1 ) lm_out <- lm_out[[1]]
    
    return( lm_out )
}
)

predict.ordlogit <- function( phi , a ) {
    K <- 1:(length(a)+1)
    a <- c( as.numeric(a) , Inf )
    p <- sapply( K , function(k) logistic(a[k]-phi) )
    na <- c(-Inf,a)
    np <- sapply( K , function(k) logistic(na[k]-phi) )
    p <- p - np
    return(p)
}
# predict.ordlogit( c(-1,0,1) , seq(-2,2,length.out=6) )

###################################
# compute WAIC from map2stan fit
# (*) needs to work with models without linear models in them

WAIC <- function( object , n=1000 , refresh=0.1 , ... ) {
    
    if ( class(object)!="map2stan" ) stop("Requires map2stan fit")
    
    if ( !is.null(attr(object,"WAIC")) ) {
        # already have it stored in object, so just return it
        return( attr(object,"WAIC") )
    }
    
    # compute linear model values at each sample
    if ( refresh > 0 ) message("Constructing posterior predictions")
    lm_vals <- link( object , n=n , refresh=refresh )
    
    # extract samples --- will need for inline parameters e.g. sigma in likelihood
    post <- extract.samples( object )
    
    n_samples <- dim( post[[1]] )[1]
    if ( n == 0 ) n <- n_samples # special flag for all samples in fit
    #if ( n_samples < n ) n <- n_samples
    
    # compute log-lik at each sample
    liks <- object@formula_parsed$likelihood
    n_lik <- length( liks )
    n_lm <- length(lm_vals)
    pD <- 0
    lppd <- 0
    flag_aggregated_binomial <- FALSE
    for ( k in 1:n_lik ) {
        outcome <- liks[[k]]$outcome
        template <- map2stan.templates[[ liks[[k]]$template ]]
        dname <- template$R_name
        pars <- liks[[k]]$pars
        # id pars as data, lm, or parameter
        pars_type <- list()
        for ( j in 1:length(pars) ) {
            pars_type[[j]] <- "UFO"
            partxt <- as.character(pars[[j]])
            if ( partxt %in% names(object@data) ) pars_type[[j]] <- "data"
            if ( partxt %in% names(post) ) pars_type[[j]] <- "parameter"
            if ( partxt %in% names(lm_vals) ) pars_type[[j]] <- "lm"
        }
        
        if ( liks[[k]]$template %in% c("Binomial") ) {
            if ( pars_type[[1]]=="data" ) flag_aggregated_binomial <- TRUE
            if ( is.numeric(pars[[1]]) )
                if ( pars[[1]]>1 ) flag_aggregated_binomial <- TRUE
        }
        if ( flag_aggregated_binomial==TRUE) {
            message("Aggregated binomial counts detected. Splitting to 0/1 outcome for WAIC calculation.")
        }
        
        #ignition
        n_obs <- liks[[k]]$N_cases
        lm_now <- list()
        
        for ( i in 1:n_obs ) {
            for ( j in 1:n_lm ) {
                # pull out samples for case i only
                ndims <- length(dim(lm_vals[[j]]))
                i_use <- min( i , dim(lm_vals[[j]])[2] )
                if ( ndims==2 )
                    lm_now[[j]] <- lm_vals[[j]][,i_use]
                if ( ndims==3 )
                    lm_now[[j]] <- lm_vals[[j]][,i_use,]
            }
            names(lm_now) <- names(lm_vals)
            
            # ready environment
            outvar <- object@data[[outcome]][i]
            # trick to get dzipois to vectorize correctly over samples
            if ( dname=="dzipois" ) outvar <- rep(outvar,n)
            e <- list( outcome=outvar )
            names(e)[1] <- outcome
            for ( j in 1:length(pars) ) {
                partxt <- as.character( pars[[j]] )
                val <- NA
                if ( pars_type[[j]]=="data" ) {
                    # fetch value for case i
                    val <- object@data[[ partxt ]][i]
                }
                if ( pars_type[[j]]=="parameter" ) {
                    # fetch n samples
                    val <- post[[ partxt ]][1:n]
                }
                if ( pars_type[[j]]=="lm" ) {
                    # fetch n values
                    ndims <- length(dim(lm_vals[[partxt]]))
                    if ( ndims==2 )
                        val <- lm_vals[[ partxt ]][ , i ]
                    if ( ndims==3 )
                        val <- lm_vals[[ partxt ]][ , i , ]
                }
                e[[ partxt ]] <- val
            }
            e1 <- new.env()
            for ( j in 1:length(e) ) {
                assign( names(e)[j] , e[[j]] , envir=e1 )
            }
            
            if ( flag_aggregated_binomial==FALSE ) {
                args_list <- list( as.symbol(outcome) , pars , log=TRUE )
                args_list <- unlist( args_list , recursive=FALSE )
                ll <- do.call( dname , args=args_list , envir=e1 )
                pD <- pD + var(ll)
                # lppd increments with log( average( likelihood ) )
                # where average is over posterior
                # but that would round to zero for sure
                # so compute on log scale with log_sum_exp
                # minus log(n) takes average over the n samples
                lppd <- lppd + log_sum_exp(ll) - log(n)
            } else {
                # aggregated binomial
                # split into 'size' 0/1 observations
                if ( is.numeric(pars[[1]]) ) 
                    size <- pars[[1]]
                if ( pars_type[[1]]=="data" ) 
                    size <- as.numeric(object@data[[as.character(pars[[1]])]][i])
                newpars <- pars
                newpars[[1]] <- 1 # bernoulli now
                # get number of 1's
                successes <- as.numeric(object@data[[as.character(outcome)]][i])
                for ( ii in 1:size ) {
                    out01 <- ifelse( ii <= successes , 1 , 0 )
                    args_list <- list( out01 , newpars , log=TRUE )
                    args_list <- unlist( args_list , recursive=FALSE )
                    ll <- do.call( dname , args=args_list , envir=e1 )
                    pD <- pD + var(ll)
                    lppd <- lppd + log_sum_exp(ll) - log(n)
                }#ii
            }#flag_aggregated_binomial
            
        }#i - cases
    }#k - linear models
    
    waic <- (-2)*( lppd - pD )
    attr(waic,"lppd") = lppd
    attr(waic,"pWAIC") = pD
    
    return(waic)
    
}




# test code
if ( FALSE ) {

library(rethinking)
data(willowtitn)
d <- willowtitn

m <- map2stan(
    alist(
        n ~ dzipois( p , lambda ),
        log(lambda) <- an + bn*percforest ,
        logit(p) ~ ap + bp*percforest ,
        c(an,ap,bn,bp) ~ dnorm(0,10)
    ),
    data=d,
    start=list(ap=0,an=1,bn=0,bp=0)
)

pred <- link( m )

d.new <- list(
    n = 1:4,
    percforest = c(10,20,30,40),
    N = 4
)
pred <- link( m , data=d.new )

######

library(rethinking)
data(UCBadmit)
d <- UCBadmit
d$male <- ifelse( d$applicant.gender=="male" , 1 , 0 )
d$dept <- as.integer(d$dept)

m <- map2stan(
    alist(
        admit ~ dbinom( applications , p ),
        logit(p) <- a + ak + b*male,
        ak[dept] ~ dnorm(0,sigma),
        a ~ dnorm(0,10),
        b ~ dnorm(0,10),
        sigma ~ dcauchy(0,2.5)
    ) ,
    data=list(admit=d$admit,applications=d$applications,male=d$male,dept=d$dept) ,
    start=list(a=0,b=0,ak=rep(0,6),sigma=1)
)

pred <- link( m )

y.mean <- apply( pred$p , 2 , mean )
y.ci <- apply( pred$p , 2 , PCI )

plot( d$admit/d$applications , ylim=c(0,1) )
lines( 1:12 , y.mean )
shade( y.ci , 1:12 )

# now with random slopes
m2 <- map2stan(
    list(
        admit ~ dbinom( applications , p ),
        logit(p) ~ a + ak + (b+bk)*male,
        c(ak,bk)|dept ~ dmvnorm(0,Sigma),
        a ~ dnorm(0,10),
        b ~ dnorm(0,10),
        Sigma ~ inv_wishart(3,diag(2))
    ) ,
    data=list(admit=d$admit,applications=d$applications,male=d$male,dept=d$dept) ,
    start=list(a=0,b=0,ak=rep(0,6),bk=rep(0,6),Sigma=diag(2))
)

pred <- link.map2stan( m2 )

y.mean <- apply( pred$p , 2 , mean )
y.ci <- apply( pred$p , 2 , PCI )

plot( d$admit/d$applications , ylim=c(0,1) )
lines( 1:12 , y.mean )
shade( y.ci , 1:12 )


#####

library(rethinking)
data(cars)

m <- map2stan(
    list(
        y ~ dnorm( mu , sigma ),
        mu ~ a + b*x,
        a ~ dnorm(0,10),
        b ~ dnorm(0,10),
        sigma ~ dcauchy(0,1)
    ) ,
    data=list(y=cars$dist,x=cars$speed) ,
    start=list(a=mean(cars$dist),b=0,sigma=10)
)

mu <- link.map2stan( m )

y.mean <- apply( mu$mu , 2 , mean )
y.ci <- apply( mu$mu , 2 , PCI )

plot( cars$dist )
lines( 1:50 , y.mean )
shade( y.ci , 1:50 )

# new data

dn <- list( y=1:50 , x=seq(0,30,length.out=50) , N=50 )

mu <- link.map2stan( m , data=dn )

y.mean <- apply( mu$mu , 2 , mean )
y.ci <- apply( mu$mu , 2 , PCI )

plot( cars$speed , cars$dist )
lines( dn$x , y.mean )
for ( p in seq(0.1,0.99,length.out=6) ) {
    y.ci <- apply( mu$mu , 2 , PCI , prob=p )
    shade( y.ci , dn$x , col=col.alpha("slateblue",0.15) )
}

}

