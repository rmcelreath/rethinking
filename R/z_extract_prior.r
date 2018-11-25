# methods for sampling from priors in map and map2stan models

setGeneric( "extract.prior",
function( fit , n=1000 , pars , ... ) {
    extract.samples(fit,n=n,pars=pars)
})

setMethod( "extract.prior", "map",
function( fit , n=1000 , pars , ... ) {

    # modified version of sample_from_prior from inside map() itself
    # this takes the formula 'language' object as input
    sample_from_prior <- function( f , n ) {
        RHS <- f[[3]]
        the_density <- as.character( RHS[[1]] )
        the_rdensity <- the_density
        substr( the_rdensity , 1 , 1 ) <- "r"
        pars <- vector(mode="list",length=length(RHS))
        pars[[1]] <- n
        for ( i in 1:(length(pars)-1) ) {
            pars[[i+1]] <- RHS[[i+1]]
            # for each argument, eval in context of data frame, in case using any constants defined there
            pars[[i+1]] <- eval( pars[[i+1]] , as.list(data) )
        }
        result <- do.call( the_rdensity , args=pars )
        return(result)
    }

    flist <- attr(fit,"formula_exploded")
    veclist <- attr(fit,"veclist")

    if ( missing(pars) ) pars <- names(fit@start)

    result <- list()
    lhs_names <- c()
    for ( i in 1:length(flist) ) {
        a_name <- flist[[i]][[2]]
        if ( length(a_name)==1 ) {
            # simple parameter
            a_name <- as.character(a_name)
        } else {
            # vector?
            if ( as.character(a_name[[1]])=="[" ) {
                # vector!
                # store just header and then use veclist in map object to get dimension later
                a_name <- as.character(a_name[[2]])
            }
        }
        lhs_names <- c( lhs_names , a_name )
    }

    for ( i in 1:length(pars) ) {

        # find formula for par
        j <- which( lhs_names == pars[i] ) # should be just one match

        if ( length(j)==0 ) {
            # vector parameter?
            pat <- strsplit( pars[i] , "[" , fixed=TRUE )[[1]][1]
            j <- grep( concat("^",pat,"\\[") , lhs_names )
        }

        # sample!
        if ( length(veclist) > 0 ) {
            if ( pars[i] %in% names(veclist) ) {
                # vector parameter, so need vector of samples
                m <- veclist[[ which( names(veclist)==pars[i] ) ]]$n
                for ( k in 1:m ) {
                    a_name <- concat( pars[i] , "[" , k , "]" )
                    result[[ a_name ]] <- sample_from_prior( flist[[j]] , n )
                }
            } else {
                result[[ pars[i] ]] <- sample_from_prior( flist[[j]] , n )
            }
        } else {
            result[[ pars[i] ]] <- sample_from_prior( flist[[j]] , n )
        }

    }#i

    # collect vector pars into an array
    name_head <- function(aname) strsplit( aname , "[" , fixed=TRUE )[[1]][1]
    name_index <- function(aname) as.numeric(regmatches( aname , regexec( "\\[(.+)\\]" , aname ) )[[1]][2])
    if ( length(veclist) > 0 ) {
        new_result <- list()
        # copy non-vector samples into new list
        for ( i in 1:length(result) ) {
            if ( !( name_head(names(result)[i]) %in% names(veclist) ) ) {
                new_result[[ names(result)[i] ]] <- result[[i]]
            }
        }#i
        for ( i in 1:length(veclist) ) {
            if ( names(veclist)[i] %in% pars ) {
                # empty matrix with parameters on cols and samples on rows
                # so n-by-m, where m in number of pars in vector and n is number of samples
                new_matrix <- matrix( 0 , ncol=veclist[[i]]$n , nrow=n )
                for ( j in 1:length(result) ) {
                    if ( name_head(names(result)[j]) == names(veclist)[i] ) {
                        the_index <- name_index( names(result)[j] )
                        new_matrix[,the_index] <- result[[j]]
                    }
                }
                new_result[[ names(veclist)[i] ]] <- new_matrix
            }
        }#i
        result <- new_result
    }

    # make sure each entry is an array, even if 1D
    for ( i in 1:length(result) ) {
        if ( is.null(dim(result[[i]])) )
            result[[i]] <- as.array(result[[i]])
    }

    return(result)
})

setMethod( "extract.prior", "map2stan",
function( fit , n=1000 , ... ) {

    # trap for ulam2018 method
    ag <- attr( fit , "generation" )
    if ( !is.null(ag) )
        if ( ag=="ulam2018" )
            return( extract_prior_ulam( fit , n=n , ... ) )

    # args must be a list of arguments to density function
    # can contain vectors
    sample_from_prior <- function( dname , args , n ) {
        the_rdensity <- dname
        substr( the_rdensity , 1 , 1 ) <- "r"
        pars <- list()
        pars$n <- n
        pars[[2]] <- args
        pars <- unlist( pars , recursive=FALSE ) #recursive=FALSE means result is a list, not vector
        result <- do.call( the_rdensity , args=pars )
        return(result)
    }

    result <- list()

    constraints <- attr(fit,"constraints")

    # sample from corresponding priors in parsed formula

    # univariate priors
    # do these first so have hyperpriors for vpriors below
    m <- length(fit@formula_parsed$prior)
    if ( m > 0 ) {
        # process in reverse so embedded parameters (lower in list) processed earlier
        for ( i in m:1 ) {
            pname <- fit@formula_parsed$prior[[i]]$par_out
            dname <- map2stan.templates[[ fit@formula_parsed$prior[[i]]$template ]]$R_name
            args <- list()
            pars_in <- fit@formula_parsed$prior[[i]]$pars_in
            for ( j in 1:length(pars_in) ) {
                if ( class(pars_in[[j]])=="numeric" ) {
                    args[[j]] <- pars_in[[j]]
                }
                if ( class(pars_in[[j]])=="name") {
                    pars_in[[j]] <- as.character(pars_in[[j]])
                }
                if ( class(pars_in[[j]])=="character" ) {
                    args[[j]] <- result[[ pars_in[[j]] ]]
                }
            }#j
            # catch rlkjcorr and insert dim K argument
            if ( dname=="dlkjcorr" ) {
                dims <- fit@stanfit@par_dims[[ pname ]]
                args <- list( dims[1] , args[[1]] )
            }
            # sample
            dims <- fit@stanfit@par_dims[[ pname ]]
            if ( length(dims)==1 ) {
                result[[pname]] <- sapply( 1:dims[1] , function(k) sample_from_prior( dname , args , n ) )
            } else
                result[[pname]] <- sample_from_prior( dname , args , n )
            # check any constraints
            if ( !is.null(constraints[[ pname ]]) ) {
                if ( constraints[[ pname ]]=="lower=0" ) {
                    # positive constraint --- hack assumes reflection at zero
                    result[[pname]] <- abs( result[[pname]] )
                }
            }
        }#i
    }

    # multivariate priors (or vector priors)
    m <- length(fit@formula_parsed$vprior)
    if ( m > 0 ) {
        for ( i in 1:m ) {
            pnames <- fit@formula_parsed$vprior[[i]]$pars_out
            dname <- map2stan.templates[[ fit@formula_parsed$vprior[[i]]$template ]]$R_name

            args <- list()
            pars_in <- fit@formula_parsed$vprior[[i]]$pars_in
            for ( j in 1:length(pars_in) ) {
                if ( class(pars_in[[j]])=="numeric" ) {
                    args[[j]] <- pars_in[[j]]
                }
                if ( class(pars_in[[j]])=="name") {
                    pars_in[[j]] <- as.character(pars_in[[j]])
                }
                if ( class(pars_in[[j]])=="character" ) {
                    args[[j]] <- result[[ pars_in[[j]] ]]
                }
            }#j
            if ( dname=="dmvnormNC" ) {
                dname <- "dmvnorm2"
                dims <- dim(args[[1]])[2] # should be length of sigma vector
                args <- list( rep(0,dims) , args[[1]] , args[[2]] )
            }

            # sample
            if ( length(pnames)==1 ) {
                # vector prior like: a[id] ~ normal( 0 , sigma_a )
                dims <- fit@stanfit@par_dims[[ pnames[[1]] ]]
                x <- sapply( 1:dims , function(k) sample_from_prior( dname , args , n ) )
                result[[ pnames[[1]] ]] <- as.array(x)
            } else {
                # multivariate like: c(a,b)[id] ~ multi_normal( 0 , Sigma )
                # sample and then shunt into correct named slots
                dims <- fit@stanfit@par_dims[[ pnames[[1]] ]] # number of units
                temp_result <- lapply( 1:dims , function(k) sample_from_prior( dname , args , n ) )
                x <- array( NA , dim=c( n , dims , length(pnames) ) )
                for ( k in 1:dims ) {
                    x[,k,] <- temp_result[[k]]
                }
                for ( k in 1:length(pnames) ) {
                    # sort columns of temp_result into named slots in result
                    result[[ pnames[[k]] ]] <- x[,,k]
                }#k
            }
        }#i
    }

    # need to get name structure from posterior samples, so can do in same order later
    post <- extract.samples(fit,n=3,...)
    new_result <- list()
    for ( i in 1:length(post) ) {
        a_name <- names(post)[i]
        if ( !is.null(result[[a_name]]) ) {
            new_result[[ a_name ]] <- result[[ a_name ]]
        }
    }#i
    result <- new_result

    # make sure each entry is an array, even if 1D
    for ( i in 1:length(result) ) {
        if ( is.null(dim(result[[i]])) )
            result[[i]] <- as.array(result[[i]])
    }

    return(result)

})

extract_prior_ulam <- function( fit , n=1000 , distribution_library=ulam_dists , ... ) {

    # args must be a list of arguments to density function
    # can contain vectors
    sample_from_prior <- function( dname , args , n ) {
        the_rdensity <- dname
        substr( the_rdensity , 1 , 1 ) <- "r"
        pars <- list()
        pars$n <- n
        pars[[2]] <- args
        pars <- unlist( pars , recursive=FALSE ) #recursive=FALSE means result is a list, not vector
        result <- do.call( the_rdensity , args=pars )
        return(result)
    }

    result <- list()

    # sample from corresponding priors in parsed formula

    # scan all symbols for parameters
    n_symbols <- length(fit@formula_parsed$symbols)
    n_pars <- 0
    pars_list <- list()
    if ( n_symbols > 0 ) {
        for ( i in 1:n_symbols ) {
            if ( fit@formula_parsed$symbols[[i]]$type=="par" ) {
                n_pars <- n_pars + 1
                par_name <- names(fit@formula_parsed$symbols)[i]
                pars_list[[ par_name ]] <- fit@formula_parsed$symbols[[i]]
            }
        }#i
    } else {
        stop( "No symbols found in formula." )
    }

    # now sample from each parameter
    if ( n_pars == 0 ) {
        stop( "No parameters found in formula." )
    }
    if ( n_pars > 0 ) {
        for ( i in 1:n_pars ) {
            
        }#i
    }

    result <- pars_list

    if ( FALSE ) {
        # need to get name structure from posterior samples, so can do in same order later
        post <- extract.samples(fit,n=3,...)
        new_result <- list()
        for ( i in 1:length(post) ) {
            a_name <- names(post)[i]
            if ( !is.null(result[[a_name]]) ) {
                new_result[[ a_name ]] <- result[[ a_name ]]
            }
        }#i
        result <- new_result

        # make sure each entry is an array, even if 1D
        for ( i in 1:length(result) ) {
            if ( is.null(dim(result[[i]])) )
                result[[i]] <- as.array(result[[i]])
        }
    }

    return(result)

}

#### TESTS #############
if ( FALSE ) {

# simple
library(rethinking)
data(cars)
flist <- alist(
    dist ~ dnorm( mu , sigma ) ,
    mu <- a+b*speed ,
    a ~ dnorm(0,0.1) , 
    b ~ dlnorm(-1,1) ,
    sigma ~ dunif(0,40)
)
fit <- map( flist , start=list(a=40,b=0.1,sigma=20) , data=cars )

prior <- extract.prior(fit,1e4)

par(mfrow=c(1,2))

mu <- link( fit , data=list(speed=0:25) , post=prior )
plot( 0:25 , apply(mu,2,mean) , ylim=c(0,max(cars$dist)) )
muci <- apply( mu , 2, PI )
shade( muci , 0:25 )
mtext("prior")

mu <- link(fit,data=list(speed=0:25))
plot( 0:25 , apply(mu,2,mean) , ylim=c(0,max(cars$dist)) )
muci <- apply( mu , 2, PI )
shade( muci , 0:25 )
mtext("posterior")

# test vectorized parameters
library(rethinking)
data(chimpanzees)

# don't want any variables with NAs
d <- list( 
    pulled_left = chimpanzees$pulled_left ,
    prosoc_left = chimpanzees$prosoc_left ,
    condition = chimpanzees$condition ,
    actor = as.integer( chimpanzees$actor ) ,
    blockid = as.integer( chimpanzees$block )
)

m <- map(
    alist(
        pulled_left ~ dbinom(1,theta),
        logit(theta) <- a + aj[actor] + bp*prosoc_left + bpc*condition*prosoc_left,
        aj[actor] ~ dnorm( 0 , 1 ),
        a ~ dnorm(0,2),
        bp ~ dnorm(0,1),
        bpc ~ dnorm(0,1)
    ) ,
    data=d )

prior <- extract.prior(m)

p <- link( m , data=data.frame(actor=1,prosoc_left=0,condition=0) , post=prior )

dens(p[,1])

prior <- extract.prior(m , pars=c("a","bp") )

prior <- extract.prior(m , pars=c("a","aj") )

# map2stan tests

m <- map2stan(
    alist(
        pulled_left ~ dbinom(1,theta),
        logit(theta) <- a + aj[actor] + bp*prosoc_left + bpc*condition*prosoc_left,
        aj[actor] ~ dnorm( 0 , sigma_actor ),
        a ~ dnorm(0,10),
        bp ~ dnorm(0,10),
        bpc ~ dnorm(0,10),
        #sigma_actor ~ dcauchy(0,1)
        sigma_actor ~ dexp(1)
    ) ,
    data=d,
    iter=1000 , chains=2 , cores=1 , sample=FALSE )

prior <- extract.prior( m , n=1e4 )
precis(prior,depth=2)
post <- extract.samples( m )
precis(post,depth=2)

par_list <- c("a","bp","bpc","sigma_actor","aj")
plot( precis(prior,2,pars=par_list),col.ci="gray")
plot( precis(post,2,pars=par_list),add=TRUE,pch=16)

m2 <- map2stan(
    alist(
        # likeliood
        pulled_left ~ dbinom(1,p),

        # linear models
        logit(p) <- A + (BP + BPC*condition)*prosoc_left,
        A <- a + a_actor[actor],
        BP <- bp + bp_actor[actor],
        BPC <- bpc + bpc_actor[actor],

        # adaptive prior - non-centered
        c(a_actor,bp_actor,bpc_actor)[actor] ~
                                dmvnormNC(sigma_actor,Rho_actor),

        # fixed priors
        c(a,bp,bpc) ~ dnorm(0,1),
        #sigma_actor ~ dcauchy(0,1),
        sigma_actor ~ dexp(1),
        Rho_actor ~ dlkjcorr(4)
    ) , data=d , iter=1000 , chains=2 , cores=1 )

prior <- extract.prior( m2 , n=1e3 )
precis(prior,depth=2)
post <- extract.samples( m2 )
precis(post,depth=2)

par_list <- c("a","bp","bpc","sigma_actor","a_actor","bp_actor","bpc_actor")
plot( precis(prior,2,pars=par_list),col.ci="gray")
plot( precis(post,2,pars=par_list),add=TRUE,pch=16)


}
