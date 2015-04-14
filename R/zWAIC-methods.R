###################################
# compute WAIC from map2stan fit
# (*) needs to work with models without linear models in them

setGeneric("WAIC",
function( object , n=1000 , refresh=0.1 , pointwise=FALSE , ... ) {
    message( concat("No WAIC method for object of class '",class(object),"'. Returning AIC instead.") )
    AIC(object)
}
)

# by default uses all samples returned by Stan; indicated by n=0
setMethod("WAIC", "map2stan",
function( object , n=0 , refresh=0.1 , pointwise=FALSE , ... ) {
    
    if ( !(class(object)%in%c("map2stan")) ) stop("Requires map2stan fit")
    
    if ( !is.null(attr(object,"WAIC")) ) {
        # already have it stored in object, so just return it
        old_waic <- attr(object,"WAIC")
        if ( length(old_waic) > 1 ) {
            if ( pointwise==FALSE ) {
                new_waic <- sum(old_waic)
                attr(new_waic,"lppd") <- sum(unlist(attr(old_waic,"lppd")))
                attr(new_waic,"pWAIC") <- sum(unlist(attr(old_waic,"pWAIC")))
                attr(new_waic,"se") <- attr(old_waic,"se")
                return( new_waic )
            } else {
                # return pointwise
                return( old_waic )
            }
        } else {
            # old waic not pointwise
            if ( pointwise==FALSE ) return( old_waic )
        }
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
    lppd_vec <- list()
    pD_vec <- list()
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
            if ( class(pars[[j]])=="call" ) {
                # language object, possible a vector of parameters in form c(...)
                if ( as.character(pars[[j]][[1]])=="c" ) pars_type[[j]] <- "parvec"
            } else {
                partxt <- as.character(pars[[j]])
                if ( partxt %in% names(object@data) ) pars_type[[j]] <- "data"
                if ( partxt %in% names(post) ) pars_type[[j]] <- "parameter"
                if ( partxt %in% names(lm_vals) ) pars_type[[j]] <- "lm"
            }
        }
        
        if ( liks[[k]]$template %in% c("Binomial") ) {
            if ( pars_type[[1]]=="data" ) flag_aggregated_binomial <- TRUE
            if ( is.numeric(pars[[1]]) )
                if ( pars[[1]]>1 ) flag_aggregated_binomial <- TRUE
        }
        if ( flag_aggregated_binomial==TRUE) {
            if ( refresh > 0 )
                message("Aggregated binomial counts detected. Splitting to 0/1 outcome for WAIC calculation.")
        }
        
        #ignition
        n_obs <- liks[[k]]$N_cases
        # vector of components so can compute std err later
        if ( flag_aggregated_binomial==FALSE ) {
            lppd_vec[[k]] <- rep(NA,n_obs)
            pD_vec[[k]] <- rep(NA,n_obs)
        } else {
            # need longer vector to hold expanded binomial
            if ( is.numeric(pars[[1]]) ) 
                size <- rep( pars[[1]] , n_obs )
            if ( pars_type[[1]]=="data" ) 
                size <- as.numeric(object@data[[as.character(pars[[1]])]])
            lppd_vec[[k]] <- rep( NA, sum(size) )
            pD_vec[[k]] <- rep( NA, sum(size) )
        }
        lm_now <- list()
        
        i_pointwise <- 1
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
                if ( pars_type[[j]]=="parvec" ) {
                    # vector of parameters
                    # for each parameter in vector, insert into environment e
                    npars <- partxt[2:length(partxt)]
                    for ( apar in 1:length(npars) ) {
                        if ( dname=="dordlogit" & !is.null(post$cutpoints) ) {
                            # get values from cutpoints
                            val <- post[[ 'cutpoints' ]][1:n,apar]
                            e[[ npars[apar] ]] <- val
                        } else {
                            val <- post[[ npars[apar] ]][1:n]
                            e[[ naprs[apar] ]] <- val
                        }
                    }
                }
                if ( pars_type[[j]]=="data" ) {
                    # fetch value for case i
                    val <- object@data[[ partxt ]][i]
                    e[[ partxt ]] <- val
                }
                if ( pars_type[[j]]=="parameter" ) {
                    # fetch n samples
                    val <- post[[ partxt ]][1:n]
                    e[[ partxt ]] <- val
                }
                if ( pars_type[[j]]=="lm" ) {
                    # fetch n values
                    ndims <- length(dim(lm_vals[[partxt]]))
                    if ( ndims==2 )
                        val <- lm_vals[[ partxt ]][ , i ]
                    if ( ndims==3 )
                        val <- lm_vals[[ partxt ]][ , i , ]
                    e[[ partxt ]] <- val
                }
            }
            e1 <- new.env()
            for ( j in 1:length(e) ) {
                assign( names(e)[j] , e[[j]] , envir=e1 )
            }
            
            if ( flag_aggregated_binomial==FALSE ) {
                args_list <- list( as.symbol(outcome) , pars , log=TRUE )
                args_list <- unlist( args_list , recursive=FALSE )
                if ( dname=="dordlogit" ) {
                    # special handling for dordlogit
                    # link() already provided probabilities for each outcome
                    # just need to extract them by indexing on outcome value
                    # 'phi' should hold matrix like phi[1:1000,1:7]
                    ll <- log( e[['phi']][ , e[[outcome]] ] )
                } else {
                    ll <- do.call( dname , args=args_list , envir=e1 )
                }
                vll <- var2(ll)
                pD <- pD + vll
                # lppd increments with log( average( likelihood ) )
                # where average is over posterior
                # but that would round to zero possibly
                # so compute on log scale with log_sum_exp
                # minus log(n) takes average over the n samples
                lpd <- log_sum_exp(ll) - log(n)
                lppd <- lppd + lpd
                # store waic for point i, so can compute stderr later
                lppd_vec[[k]][i] <- lpd
                pD_vec[[k]][i] <- vll
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
                    vll <- var2(ll)
                    pD <- pD + vll
                    lpd <- log_sum_exp(ll) - log(n)
                    lppd <- lppd + lpd
                    lppd_vec[[k]][i_pointwise] <- lpd
                    pD_vec[[k]][i_pointwise] <- vll
                    i_pointwise <- i_pointwise + 1
                }#ii
            }#flag_aggregated_binomial
            
        }#i - cases
    }#k - linear models
    
    waic_vec <- (-2)*( unlist(lppd_vec) - unlist(pD_vec) )
    if ( pointwise==TRUE ) {
        waic <- unlist(waic_vec)
        lppd <- unlist(lppd_vec)
        pD <- unlist(pD_vec)
    } else {
        waic <- (-2)*( lppd - pD )
    }
    attr(waic,"lppd") = lppd
    attr(waic,"pWAIC") = pD
    
    n_tot <- length(waic_vec)
    attr(waic,"se") = try(sqrt( n_tot*var2(waic_vec) ))
    
    return(waic)
    
}
)

setMethod("WAIC", "map",
function( object , n=1000 , refresh=0.1 , pointwise=FALSE , ... ) {
    
    if ( !(class(object)%in%c("map")) ) stop("Requires map fit")
    
    if ( !is.null(attr(object,"WAIC")) ) {
        # already have it stored in object, so just return it
        old_waic <- attr(object,"WAIC")
        if ( pointwise==TRUE ) {
            if ( length(old_waic)>1 ) {
                # already have pointwise version, so return it
                return( old_waic )
            } else {
                # old one is not pointwise, so need to recalcuate
                ### DO NOTHING HERE
            }
        } else {
            if ( length(old_waic)>1 ) {
                # old one was pointwise, so total it up
                new_waic <- sum(old_waic)
                attr(new_waic,"lppd") <- sum(attr(old_waic,"lppd"))
                attr(new_waic,"pWAIC") <- sum(attr(old_waic,"pWAIC"))
                attr(new_waic,"se") <- attr(old_waic,"se")
            }
        }
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
    lppd_vec <- rep(NA,n_cases)
    pD_vec <- rep(NA,n_cases)
    for ( i in 1:n_cases ) {# for each case
        
        vll <- var2(s[,i])
        pD <- pD + vll
        lpd <- log_sum_exp(s[,i]) - log(n_samples)
        lppd <- lppd + lpd
        lppd_vec[i] <- lpd
        pD_vec[i] <- vll
    }#i - cases
    
    waic_vec <- (-2)*( lppd_vec - pD_vec )
    if ( pointwise==TRUE ) {
        # return decomposed as WAIC for each observation i --- can sum to get total WAIC
        # this is useful for computing standard errors of differences with compare()
        waic <- (-2)*( lppd_vec - pD_vec )
        lppd <- lppd_vec
        pD <- pD_vec
    } else {
        waic <- (-2)*( lppd - pD )
    }
    attr(waic,"lppd") = lppd
    attr(waic,"pWAIC") = pD
    attr(waic,"se") = try(sqrt( n_cases*var2(waic_vec) ))
    
    return(waic)
    
}
)

# lm
setMethod("WAIC", "lm",
function( object , n=1000 , refresh=0.1 , pointwise=FALSE , ... ) {
    
    if ( refresh > 0 ) message("Constructing posterior predictions")
    
    # extract samples --- will need for inline parameters e.g. sigma in likelihood
    post <- extract.samples( object , n=n )
    # get sigma
    sigma <- sd(resid(object))
    
    # compute log-lik at each sample
    # can hijack predict.lm() to do this by inserting samples into object$coefficients
    m_temp <- object
    out_var <- object$model[[1]] # should be outcome variable
    n_coefs <- length(coef(object))
    s <- sapply( 1:n , 
        function(i) {
            for ( j in 1:n_coefs ) m_temp$coefficients[[j]] <- post[i,j]             
            mu <- as.numeric(predict(m_temp)) # as.numeric strips names
            dnorm( out_var , mu , sigma , log=TRUE )
        } )
    
    pD <- 0
    lppd <- 0
    n_cases <- length(out_var)
    n_samples <- n
    lppd_vec <- rep(NA,n_cases)
    pD_vec <- rep(NA,n_cases)
    for ( i in 1:n_cases ) {# for each case
        
        vll <- var2(s[i,])
        pD <- pD + vll
        lpd <- log_sum_exp(s[i,]) - log(n_samples)
        lppd <- lppd + lpd
        lppd_vec[i] <- lpd
        pD_vec[i] <- vll
    }#i - cases
    
    waic_vec <- (-2)*( lppd_vec - pD_vec )
    if ( pointwise==TRUE ) {
        # return decomposed as WAIC for each observation i --- can sum to get total WAIC
        # this is useful for computing standard errors of differences with compare()
        waic <- (-2)*( lppd_vec - pD_vec )
        lppd <- lppd_vec
        pD <- pD_vec
    } else {
        waic <- (-2)*( lppd - pD )
    }
    attr(waic,"lppd") = lppd
    attr(waic,"pWAIC") = pD
    attr(waic,"se") = try(sqrt( n_cases*var2(waic_vec) ))
    
    return(waic)
} )

