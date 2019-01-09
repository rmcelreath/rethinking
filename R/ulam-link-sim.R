
# new link function
link_ulam <- function( fit , data , post , simplify=TRUE , symbols , ... ) {

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

    for ( j in 1:n_links ) {

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
                if ( j > 1 ) {
                    # insert previous linear model results, in case used in later ones
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

    n_links <- length(symbols)

    if ( simplify==TRUE && n_links==1 ) out <- out[[1]]

    return(out)
}
setMethod( "link" , "ulam" , function(fit ,data ,...) link_ulam(fit,data,...) )

# sim method assumes outcome var is first var on left in formula line 1
# can use variable argument to specify any other variable
sim_ulam <- function( fit , data , post , variable , n=1000 , replace=list() , ... ) {
    
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
    for ( s in 1:n ) {
        # build environment
        #ndims <- length(dim(pred[[1]]))
        #if ( ndims==2 )
        #    elm <- pred[[1]][s,] # extract s-th sample case calculations
        #if ( ndims==3 )
        #    elm <- pred[[1]][s,,]
        #elm <- list(elm)
        #names(elm)[1] <- names(pred)[1]
        
        elm <- list()
        for ( j in 1:length(pred) ) {
            ndims <- length(dim(pred[[j]]))
            if ( ndims==2 )
                elm[[j]] <- pred[[j]][s,]
            if ( ndims==3 )
                elm[[j]] <- pred[[j]][s,,]
        }
        names(elm) <- names(pred)
        
        e <- list( as.list(data) , as.list(post2[s,]) , as.list(elm) )
        e <- unlist( e , recursive=FALSE )
        # evaluate
        if ( flik=="dordlogit" ) {
            n_outcome_vals <- dim( pred[[1]] )[3]
            probs <- pred[[1]][s,,]
            sim_out[s,] <- sapply( 
                1:n_cases ,
                function(i)
                    sample( 1:n_outcome_vals , size=1 , replace=TRUE , prob=probs[i,] )
            )
        } else {
            sim_out[s,] <- eval(parse(text=xeval),envir=e)
        }
    }
    
    return(sim_out)
}
setMethod( "sim" , "ulam" , function(fit,...) sim_ulam(fit,...) )


