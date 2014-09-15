###################################
# compute WAIC from map2stan fit
# (*) needs to work with models without linear models in them

setGeneric("WAIC",
function( object , n=1000 , refresh=0.1 , ... ) {
    message( concat("No WAIC method for object of class '",class(object),"'. Returning AIC instead.") )
    AIC(object)
}
)

setMethod("WAIC", "map2stan",
function( object , n=1000 , refresh=0.1 , pointwise=FALSE , ... ) {
    
    if ( !(class(object)%in%c("map2stan")) ) stop("Requires map2stan fit")
    
    if ( !is.null(attr(object,"WAIC")) ) {
        # already have it stored in object, so just return it
        old_waic <- attr(object,"WAIC")
        if ( length(old_waic) > 1 & pointwise==FALSE ) {
            old_waic <- (-2)*sum(old_waic)
            return( old_waic )
        }
        if ( length(old_waic) > 1 & pointwise==TRUE ) {
            return( old_waic )
        } # else continue onward to recompute
    }
    
    # compute linear model values at each sample
    if ( refresh > 0 ) message("Constructing posterior predictions")
    lm_vals <- link( object , n=n , refresh=refresh , flatten=FALSE )
    
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
    waic_vec <- list()
    flag_aggregated_binomial <- FALSE
    for ( k in 1:n_lik ) {
        outcome <- liks[[k]]$outcome
        
        # check for predictor imputation definition
        #   if so, skip this likelihood, because not really an outcome to predict
        if ( outcome %in% names(object@formula_parsed$impute_bank) ) {
            # message( concat("skipping ",outcome) )
            next
        }
        
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
        waic_vec[[k]] <- rep(NA,n_obs) # vector of components so can compute std err later
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
                vll <- var(ll)
                pD <- pD + vll
                # lppd increments with log( average( likelihood ) )
                # where average is over posterior
                # but that would round to zero possibly
                # so compute on log scale with log_sum_exp
                # minus log(n) takes average over the n samples
                lpd <- log_sum_exp(ll) - log(n)
                lppd <- lppd + lpd
                # store waic for point i, so can compute stderr later
                waic_vec[[k]][i] <- lpd - vll
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
                    vll <- var(ll)
                    pD <- pD + vll
                    lpd <- log_sum_exp(ll) - log(n)
                    lppd <- lppd + lpd
                    waic_vec[[k]][i] <- lpd - vll
                }#ii
            }#flag_aggregated_binomial
            
        }#i - cases
    }#k - linear models
    
    if ( pointwise==TRUE ) {
        waic <- unlist(waic_vec)
    } else {
        waic <- (-2)*( lppd - pD )
    }
    attr(waic,"lppd") = lppd
    attr(waic,"pWAIC") = pD
    
    waic_vec <- unlist(waic_vec)
    n_tot <- length(waic_vec)
    attr(waic,"se") = try(sqrt( n_tot*var(waic_vec) ))
    
    return(waic)
    
}
)

setMethod("WAIC", "map",
function( object , n=1000 , refresh=0.1 , pointwise=FALSE , ... ) {
    
    if ( !(class(object)%in%c("map")) ) stop("Requires map fit")
    
    if ( !is.null(attr(object,"WAIC")) ) {
        # already have it stored in object, so just return it
        return( attr(object,"WAIC") )
    }
    
    # compute linear model values at each sample
    if ( refresh > 0 ) message("Constructing posterior predictions")
    #lm_vals <- link( object , n=n , refresh=refresh , flatten=FALSE )
    
    # extract samples --- will need for inline parameters e.g. sigma in likelihood
    post <- extract.samples( object , n=n )
    
    # compute log-lik at each sample
    #lik <- flist_untag(object@formula)[[1]]
    s <- sim(object,post=post,ll=TRUE,refresh=refresh)
    pD <- 0
    lppd <- 0
    flag_aggregated_binomial <- FALSE
    n_cases <- dim(s)[2]
    n_samples <- dim(s)[1]
    waic_vec <- rep(NA,n_cases)
    for ( i in 1:n_cases ) {# for each case
        
        vll <- var(s[,i])
        pD <- pD + vll
        lpd <- log_sum_exp(s[,i]) - log(n_samples)
        lppd <- lppd + lpd
        waic_vec[i] <- lpd - vll
        
    }#i - cases
    
    if ( pointwise==TRUE ) {
        # return decomposed as WAIC for each observation i --- can sum to get total WAIC
        # this is useful for computing standard errors of differences with compare()
        waic <- waic_vec
    } else {
        waic <- (-2)*( lppd - pD )
    }
    attr(waic,"lppd") = lppd
    attr(waic,"pWAIC") = pD
    attr(waic,"se") = try(sqrt( n_cases*var(waic_vec) ))
    
    return(waic)
    
}
)
