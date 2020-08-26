# sim -- simulates observations from fit map/quap and map2stan and ulam models

setGeneric("sim",
function( fit , data , n=1000 , ... ) {
    predict(fit)
}
)

# function at heart of all sim methods
sim_core <- function( fit , data , post , vars , n , refresh=0 , replace=list() , debug=FALSE , ll=FALSE , ... ) {

    sim_vars <- list()

    for ( var in vars ) {

        # find a distribtional assignment with var on left
        f <- fit@formula[[1]]
        for ( i in 1:length(fit@formula) ) {
            f <- fit@formula[[i]]
            left <- as.character( f[[2]] )
            if ( left==var ) {
                if ( debug==TRUE ) print(f)
                break
            }
        }

        lik <- f
        outcome <- var
        flik <- as.character(lik[[3]][[1]])
        # check whether we must convert from Stan-name distribution to R-name distribution
        first_char <- substr( flik , 1 , 1 )
        if ( first_char != "d" | flik=="dirichlet" ) {
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
        if ( ll==FALSE ) {
            substr( rlik , 1 , 1 ) <- "r"
            # hook for gampois, which needs rgampois2
            if ( rlik=="rgampois" ) rlik <- "rgampois2"
        }
        
        # check for aggregated binomial, but only when ll==TRUE
        aggregated_binomial <- FALSE
        size_var_is_data <- FALSE
        if ( flik=="dbinom" & ll==TRUE ) {
            ftemp <- flist_untag(fit@formula)[[1]]
            size_sym <- ftemp[[3]][[2]]
            if ( class(size_sym)=="name" ) {
                aggregated_binomial <- TRUE
                size_var_is_data <- TRUE
            }
            if ( class(size_sym)=="numeric" ) {
                if ( size_sym > 1 ) aggregated_binomial <- TRUE
            }
        }

        # disabling aggregated binom code for now
        aggregated_binomial <- FALSE
        
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
        if ( n_cases==0 ) {
            # no outcome in custom data?
            # get number of cases from first variable in data
            n_cases <- length(data[[1]])
        }
        if ( ll==FALSE ) {
            xeval <- paste( rlik , "(" , n_cases , "," , pars , ")" , collapse="" )
        } else {
            use_outcome <- outcome
            if ( aggregated_binomial==TRUE ) use_outcome <- 'ones__'
            xeval <- paste( rlik , "(" , use_outcome , "," , pars , ",log=TRUE )" , collapse="" )
        }
        
        # simulate outcomes
        sim_out <- matrix(NA,nrow=n,ncol=n_cases)

        # get linear model values from link
        # don't flatten result, so we end up with a named list, even if only one element
        # need to insert already-simulated vars here - should have happened at bottom of this loop
        pred <- link( fit , data=data , n=n , post=post , flatten=FALSE , refresh=refresh , replace=replace , debug=debug , ... )

        for ( s in 1:n ) {
            # build environment
            # extract s-th sample case calculations
            elm <- list()
            for ( j in 1:length(pred) ) {
                elm[[j]] <- pred[[j]][s,]
            }
            names(elm) <- names(pred)
            init <- list() # holds one row of samples across all params
            for ( j in 1:length(post) ) {
                par_name <- names(post)[ j ]
                dims <- dim( post[[par_name]] )
                # scalar
                if ( is.null(dims) ) init[[par_name]] <- post[[par_name]][s]
                # 1d vector
                if ( length(dims)==1 ) init[[par_name]] <- post[[par_name]][s]
                # vector
                if ( length(dims)==2 ) init[[par_name]] <- post[[par_name]][s,]
                # matrix
                if ( length(dims)==3 ) init[[par_name]] <- post[[par_name]][s,,]
            }#j

            e <- list( as.list(data) , as.list(init) , as.list(elm) )
            e <- unlist( e , recursive=FALSE )
            # include previously simulated variables in environment so they are accessible
            # need only s-th simulation right now, and should replace existing symbols in e (in case variables were included in data list)
            if ( length(sim_vars) > 0 ) {
                for ( j in 1:length(sim_vars) ) {
                    e[[ names(sim_vars)[j] ]] <- sim_vars[[j]][s,]
                }#j
            }

            # evaluate
            sim_out[s,] <- eval( parse( text=xeval ) , envir=e )
        }#s
        if ( debug==TRUE ) print(str(e))
        
        # check for aggregated binomial outcome
        if ( aggregated_binomial==TRUE ) {
            # aggregated binomial with data for 'size'
            # need to split binomial counts in outcome into series of 0/1 outcomes
            # (1) sum 'size' variable in order to get number of new cases
            if ( size_var_is_data==TRUE ) {
                size_var <- data[[ as.character(ftemp[[3]][[2]]) ]]
            } else {
                # size is constant numeric, but > 1
                size_var <- rep( as.numeric(ftemp[[3]][[2]]) , n_cases )
            }
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
                if ( num_ones > 0 ) {
                    for ( j in 1:num_ones ) {
                        sim_out_new[,current_newcase] <- ll1
                        current_newcase <- current_newcase + 1
                    }
                }
                if ( num_zeros > 0 ) {
                    for ( j in 1:num_zeros ) {
                        sim_out_new[,current_newcase] <- ll0
                        current_newcase <- current_newcase + 1
                    }
                }
            }#i
            sim_out <- sim_out_new
        }

        sim_vars[[ var ]] <- sim_out
        data[[ var ]] <- sim_out # make available to subsequent vars

    }#var

    return(sim_vars)

}

setMethod("sim", "map",
function( fit , data , n=1000 , post , vars , ll=FALSE , refresh=0 , replace=list() , debug=FALSE , ... ) {
    # when ll=FALSE, simulates sampling, one sample for each sample in posterior
    # when ll=TRUE, computes loglik of each observation, for each sample in posterior

    # vars argument optional
    # when missing, simulations only first observed variable in formula
    # when present, should be a character vector naming variables to simulate, in order to simulate them
    # order matters, because some may be needed for the others to simulate right
    
    ########################################
    # check arguments
    if ( class(fit)!="map" ) stop("Requires map/quap fit")
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
    
    if ( missing(post) ) {
        post <- extract.samples(fit,n=n)
    } else {
        n <- dim(post[[1]])[1]
        if ( is.null(n) ) n <- length(post[[1]])
    }
    
    # variables to sim
    if ( missing(vars) ) {
        # extract likelihood, assuming it is first element of formula
        lik <- flist_untag(fit@formula)[[1]]
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
)

setMethod("sim", "map2stan",
function( fit , data , n=1000 , post , refresh=0.1 , replace=list() , ... ) {
    
    # trap for ulam2018 method
    ag <- attr( fit , "generation" )
    if ( !is.null(ag) )
        if ( ag=="ulam2018" )
            return( sim_ulam( fit , data=data , post=post , n=n , ... ) )

    ########################################
    # check arguments
    if ( missing(data) ) {
        data <- fit@data
    } else {
        # make sure all variables same length
        # weird vectorization errors otherwise
        #data <- as.data.frame(data)
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
    pred <- link( fit , data=data , n=n , post=post , flatten=FALSE , refresh=refresh , replace=replace , ... )
    
    # extract likelihood, assuming it is first element of formula
    lik <- flist_untag(fit@formula)[[1]]
    # discover outcome
    outcome <- as.character(lik[[2]])
    # discover likelihood function
    flik <- as.character(lik[[3]][[1]])

    # check whether we must convert from Stan-name distribution to R-name distribution
    first_char <- substr( flik , 1 , 1 )
    if ( first_char != "d" ) {
        # loop over templates and get matching dfoo name
        for ( ii in 1:length(map2stan.templates) ) {
            aStanName <- map2stan.templates[[ii]]$stan_name
            if ( aStanName==flik ) {
                flik <- map2stan.templates[[ii]]$R_name
                break
            }
        }#ii
    }

    # get simulation partner function
    rlik <- flik
    substr( rlik , 1 , 1 ) <- "r"
    # hook for gampois, which needs rgampois2
    if ( rlik=="rgampois" ) rlik <- "rgampois2"

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
)

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
