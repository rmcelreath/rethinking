
# new link function
link_ulam <- function( fit , data , post , flatten=TRUE , symbols , ... ) {

    nest_commas <- function( f , symbols ) {
        if ( class( f )=="call" || class( f )=="(" ) {
            if ( as.character( f[[1]] )=="[" ) {
                # bracket call, so maybe insert comma (TRUE)
                if ( as.character( f[[2]] ) %in% symbols ) {
                    the_call <- as.list( f )
                    f <- as.call( unlist( list( the_call[1:2] , quote(TRUE) , the_call[3:length(the_call)] ) , recursive=FALSE ) )
                }
            } else {
                # need to drill down
                for ( i in 1:length(f) )
                    f[[i]] <- nest_commas( f[[i]] , symbols )
            }
        }
        return( f )
    }

    use_orig_data <- FALSE
    if ( missing(data) ) {
        data <- fit@data
        use_orig_data <- TRUE
    }
    if ( missing(post) ) post <- extract.samples(fit@stanfit)

    # how many link functions?
    if ( !is.null(fit@formula_parsed$link_funcs) ) {
        n_links <- length( fit@formula_parsed$link_funcs )
    } else {
        stop( "No link functions found. Maybe there were no linear models in the formula list?" )
    }

    out <- list()

    if ( missing(symbols) ) {
        # use them all
        symbols <- names( fit@formula_parsed$link_funcs )
    }

    # go backwards, so links further down in formula can be insterted above
    for ( j in n_links:1 ) {

        if ( !(names( fit@formula_parsed$link_funcs )[j] %in% symbols) ) next

        # get number of cases for this symbol
        if ( use_orig_data==TRUE )
            n_cases <- fit@formula_parsed$link_funcs[[ j ]]$N
        else
            # guess from length of data
            n_cases <- length(data[[1]])

        # check whether this function contains symbols that need some dimensions added for samples
        symbols_so_far <- names(out)
        f <- fit@formula_parsed$link_funcs[[ j ]]$func
        for ( jj in symbols_so_far ) {
            # more than one dim? if so, then may need extra dim when referenced
            if ( length(dim( out[[jj]] ))>1 ) {
                f <- nest_commas( f , jj )
            }
        }

        l <- sapply( 1:n_cases , function(i) {
                the_env <- post
                #the_env <- unlist( list( the_env , data[i,,drop=FALSE] ) , recursive=FALSE )
                the_env <- unlist( list( the_env , data , i=i ) , recursive=FALSE )
                if ( j < n_links ) {
                    # insert previous linear model results, in case used in later ones
                    if ( length(out)>0 )
                        for ( k in 1:length(out) ) {
                            lv_name <- names(out)[k]
                            #the_env[[ lv_name ]] <- out[[k]][ , i ] # all samples, i-th case
                            the_env[[ lv_name ]] <- out[[k]]
                        }
                }
                eval( f , envir=the_env ) 
            })

        # an inverse link to apply?
        if ( !is.null( fit@formula_parsed$link_funcs[[ j ]]$inv_link ) ) {
            l <- do.call( fit@formula_parsed$link_funcs[[ j ]]$inv_link , list(l) )
        }

        out[[ names( fit@formula_parsed$link_funcs )[j] ]] <- l

    }#j

    # reverse order of elements in out
    outn <- list()
    for ( j in length(out):1 ) outn[[ names(out)[j] ]] <- out[[ j ]]
    out <- outn

    n_links <- length(symbols)

    if ( flatten==TRUE && n_links==1 ) out <- out[[1]]

    return(out)
}
setMethod( "link" , "ulam" , function(fit ,data ,...) link_ulam(fit,data,...) )

# sim method assumes outcome var is first var on left in formula line 1
# can use variable argument to specify any other variable
sim_ulam <- function( fit , data , post , vars , variable , n=1000 , replace=list() , ... ) {
    
    ########################################
    # check arguments
    if ( missing(data) ) {
        data <- fit@data
    }

    if ( n==0 ) {
        n <- stan_total_samples(fit@stanfit)
    } else {
        tot_samples <- stan_total_samples(fit@stanfit)
        n <- min(n,tot_samples)
    }
    
    if ( missing(post) ) 
        post <- extract.samples(fit,n=n)

    # get linear model values from link
    # use our posterior samples, so later parameters have right correlation structure with link values
    # don't flatten result, so we end up with a named list, even if only one element
    pred <- link_ulam( fit , data=data , post=post , simplify=FALSE , ... )
    
    # extract likelihood, assuming it is first element of formula
    if ( missing(variable) ) {
        lik <- fit@formula[[1]]
        # discover outcome
        outcome <- as.character(lik[[2]])
        # discover likelihood function
        flik <- as.character(lik[[3]][[1]])
    } else {
        # named outcome, so find it
        for ( i in 1:length(fit@formula) ) {
            outcome <- as.character( fit@formula[[i]][[2]] )
            if ( outcome==variable ) {
                lik <- fit@formula[[i]]
                flik <- as.character( lik[[3]][[1]] )
                break
            }
        }#i
    }

    # check whether we must convert from Stan-name distribution to R-name distribution
    first_char <- substr( flik , 1 , 1 )
    if ( first_char != "d" ) {
        # loop over templates and get matching dfoo name
        for ( ii in 1:length( ulam_dists ) ) {
            aStanName <- ulam_dists[[ii]]$Stan_name
            if ( aStanName==flik ) {
                flik <- ulam_dists[[ii]]$R_name
                break
            }
        }#ii
    }

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
    post2 <- as.data.frame(post)
    cuts_name <- "cutpoints"
    for ( s in 1:n ) {
        
        # handle ordered logit as special case
        # it is very fragile
        if ( flik=="dordlogit" ) {
            n_outcome_vals <- dim( post[[ cuts_name ]] )[2] + 1
            phi <- pred[[1]][s,]
            #probs <- pordlogit( 1:n_outcome_vals , phi , post[[ cuts_name ]][s,] )
            #sim_out[s,] <- sapply( 
            #    1:n_cases ,
            #    function(i)
            #        sample( 1:n_outcome_vals , size=1 , replace=TRUE , prob=probs[i,] )
            #)
            sim_out[s,] <- rordlogit( length(phi) , phi , post[[ cuts_name ]][s,] )
        } else {
            # build environment
            elm <- list()
            for ( j in 1:length(pred) ) {
                ndims <- length(dim(pred[[j]]))
                if ( ndims==2 )
                    elm[[j]] <- pred[[j]][s,]
                if ( ndims==3 )
                    elm[[j]] <- pred[[j]][s,,]
            }
            names(elm) <- names(pred)
            # evaluate
            e <- list( as.list(data) , as.list(post2[s,]) , as.list(elm) )
            e <- unlist( e , recursive=FALSE )
            sim_out[s,] <- eval(parse(text=xeval),envir=e)
        }
    }
    
    return(sim_out)
}

sim_ulam_new <- function( fit , data , post , vars , variable , n=1000 , replace=list() , debug=FALSE , ll=FALSE , refresh=0 , ... ) {
    
    ########################################
    # check arguments
    if ( missing(data) ) {
        data <- fit@data
    } else {
        # make sure it's a list, otherwise can't hold sample matrices from sims
        data <- as.list(data)
        # check for vars to sim that are not in data list
        if ( !missing(vars) ) {
            for ( i in 1:length(vars) ) {
                if ( !(vars[i] %in% names(data)) ) {
                    # insert dummy for var - will get simulated over later
                    data[[ vars[i] ]] <- rep( 0 , length(data[[1]]) )
                }
            }#i
        }
    }

    if ( n==0 ) {
        n <- stan_total_samples(fit@stanfit)
    } else {
        tot_samples <- stan_total_samples(fit@stanfit)
        n <- min(n,tot_samples)
    }
    
    if ( missing(post) ) {
        post <- extract.samples(fit,n=n)
    } else {
        n <- dim(post[[1]])[1]
        if ( is.null(n) ) n <- length(post[[1]])
    }

    # variables to sim
    if ( missing(vars) ) {
        # extract likelihood, assuming it is first element of formula
        #lik <- flist_untag(fit@formula)[[1]]
        lik <- fit@formula[[1]]
        # discover outcome
        vars <- as.character(lik[[2]])
    } else {
        # vars listed
        if ( debug==TRUE ) print(vars)
    }

    # loop over vars
    sim_vars <- sim_core( fit=fit , data=data , post=post , vars=vars , n=n , refresh=refresh , replace=replace , debug=debug , ll=ll , ... )
    
    # result
    if ( length(sim_vars)==1 ) sim_vars <- sim_vars[[1]]
    return( sim_vars )

}

setMethod( "sim" , "ulam" , function(fit,data,...) sim_ulam_new(fit,data,...) )


