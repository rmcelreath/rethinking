# MAP - maximum a posteriori estimation

# to-do:

# utility function to convert alist() construction with <- tagged linear model to regular ~ variety
# do some error checking here
flist_untag <- function(flist) {
    for ( i in 1:length(flist) ) {
        if ( !is.null(names(flist)) ) {
            if ( names(flist)[i]!="" ) {
                warning( concat("Named entry '",names(flist)[i],"' detected. Make sure you didn't use '=' where you meant '~' or '<-'.") )
            }
        }#!is.null
        if ( class(flist[[i]])=="<-" ) {
            # tagged formula, so convert to ~ formula expression
            flist[[i]][[1]] <- as.name("~")
        }
        # eval required to convert from class language to class formula
        flist[[i]] <- eval(flist[[i]])
    }
    as.list(flist)
}

# main event
map <- function( flist , data , start , method="BFGS" , hessian=TRUE , debug=FALSE , verbose=FALSE , ... ) {
    
    ########################################
    # check arguments
    if ( missing(flist) ) stop( "Formula required." )
    if ( class(flist) != "list" ) {
        if ( class(flist)=="formula" ) {
            flist <- list(flist)
        } else {
            stop( "Formula or list of formulas required." )
        }
    }
    if ( missing(start) ) start <- list()
    if ( missing(data) ) stop( "'data' required." )
    if ( !( class(data) %in% c("list","data.frame") ) ) {
        stop( "'data' must be of class list or data.frame." )
    }
    
    flist.orig <- flist
    flist <- flist_untag(flist)
    
    ########################################
    # check for common issues
    # (1) NAs in data vectors
    # tricky, because only want to check used variables
    
    # (2) type mismatch for outcome variable
    
    ########################################
    # private functions
    
    # substitution list - NOT YET IMPLEMENTED
    # will allow common (Stan-like) names for distributions
    density_sub_list <- list(
        normal = 'dnorm',
        binomial = 'dbinom',
        poisson = 'dpois',
        gamma = 'dgamma2', # mu,scale version
        betabinomial = 'dbetabinom', #in emdbook package
        gammapoisson = 'dgampois',
        cauchy = 'dcauchy',
        uniform = 'dunif',
        laplace = 'dlaplace'
    )
    
    idx_marker_string <- "___"
    
    # function to sample from prior specified with density function
    sample_from_prior <- function( f ) {
        RHS <- f[[3]]
        the_density <- as.character( RHS[[1]] )
        the_rdensity <- the_density
        substr( the_rdensity , 1 , 1 ) <- "r"
        pars <- vector(mode="list",length=length(RHS))
        pars[[1]] <- 1
        for ( i in 1:(length(pars)-1) ) {
            pars[[i+1]] <- RHS[[i+1]]
            # for each argument, eval in context of data frame, in case using any constants defined there
            pars[[i+1]] <- eval( pars[[i+1]] , as.list(data) )
        }
        result <- do.call( the_rdensity , args=pars )
        return(result)
    }
    
    # these functions are the engine that does the optimization
    # trick is to build R code for objective function that is passed to parser
    # it's slow, but easier to debug

    dparser <- function( flist , e ) {
    # merge this into make_minuslogl?
        r <- sapply( flist , function(i) sum(eval(parse(text=i),envir=e)) )
        -sum(r)
    }
    
    # this function is needed to make optim work with vector parameters
    pars_to_vectors <- function( pars_orig , veclist ) {
        pars_new <- list()
        # first copy all scalar parameters
        for ( i in 1:length(pars_orig) ) {
            x <- strsplit( names(pars_orig)[i] , idx_marker_string )[[1]]
            if ( length(x)==1 ) {
                pars_new[[ x ]] <- pars_orig[[i]]
            }
        }
        # now compose vectors
        for ( i in 1:length(veclist) ) {
            newvec <- rep( 0 , veclist[[i]]$n )
            for ( j in 1:veclist[[i]]$n ) {
                name_orig <- paste( veclist[[i]]$name , idx_marker_string , j , sep="" , concat="" )
                newvec[j] <- pars_orig[[ name_orig ]]
            }
            pars_new[[ veclist[[i]]$name ]] <- newvec
        }#i
        return(pars_new)
    }
    
    make_minuslogl <- function( pars , flist , data , veclist ) {
        # check for pars with name that ends in _._n
        # convert them to vector
        if ( length(veclist)>0 ) pars <- pars_to_vectors( pars , veclist )
        e <- list( as.list(data) , as.list(pars) )
        e <- unlist( e , recursive=FALSE )
        dparser(flist,e)
    }
    
    link.names <- c("log","logit")
    invlink.names <- c("exp","logistic")
    formula2text <- function( f ) {
    # must keep named arguments with names, 
    # so user can pass arguments out of order to density function
        RHS <- f[[3]]
        LHS <- f[[2]]
        flag_monad_linear_model <- FALSE
        if ( length(RHS)==1 ) {
            if ( class(RHS)=="numeric" | class(RHS)=="name" )
                flag_monad_linear_model <- TRUE
            fname <- ""
        } else {
            fname <- as.character( RHS[[1]] )
        }
        if ( fname=="[" | fname=="+" | fname=="*" | fname=="-" | fname=="/" | fname=="%*%" | fname %in% invlink.names | flag_monad_linear_model==TRUE ) {
            # linear model formula with no density (but maybe invlink) function
            # return a list with parameter name in [[1]] and text of RHS in [[2]]
            thetext <- list( as.character(LHS) , paste( deparse(RHS) , collapse=" " ) )
        } else {
            # likelihood or prior formula
            n_args <- length(RHS)
            args_list <- as.list(RHS)
            # check for LHS with brackets []
            if ( class(LHS)=="call" ) {
                if ( as.character(LHS[[1]])=="[" ) {
                    ival <- suppressWarnings( as.numeric(as.character(LHS[[3]])) )
                    if ( is.na(ival) ) {
                        # [var] form, so just delete the bracket so prior properly vectorizes over vector of parameters, not over entire data table
                        LHS <- as.character(LHS[[2]])
                    } else {
                        # [num] form, so preserve index number, because prior is specific to this index
                        LHS <- deparse(LHS)
                    }
                }
            }
            args_list[[1]] <- LHS
            args_list[[n_args+1]] <- "TRUE"
            args_names <- names(RHS)[-1]
            if ( is.null(args_names) ) args_names <- rep("",n_args-1)
            args_names <- c( "x" , args_names , "log" )
            for ( i in 1:length(args_names) ) {
                if ( args_names[i]=="" ) {
                    # no argument name
                    args_names[i] <- args_list[i]
                } else {
                    # prepend argument name
                    args_names[i] <- paste( args_names[i] , args_list[i] , sep="=" )
                }
            }
            args_text <- paste( args_names , collapse=" , " )
            thetext <- paste( RHS[[1]] , "(" , args_text , ")" , sep="" )
        }
        return( thetext )
    }
    
    ##############
    # grep function for search-replace of linear model symbols
    # trick is that symbol needs to be preceded by [(=*+ ] so grep doesn't replace copies embedded inside other symbols
    # e.g. don't want to expand the "p" inside "ample"
    #
    # target : string to search for (usually parameter name like "mu")
    # replacement : what to replace with (usually a linear model)
    # x : where to search, usually a formula as character
    # add.par : whether to enclose replacement in parentheses
    
    # needed , in case of vector or list construction in likelihood
    mygrep <- function( target , replacement , x , add.par=TRUE ) {
        wild <- "[()=*+, ]"
        pattern <- paste( wild , target , wild , sep="" , collapse="" )
        m <- regexpr( pattern , x )
        if ( m==-1 ) return( x )
        s <- regmatches( x=x , m=m )
        
        if ( add.par==TRUE ) replacement <- paste( "(" , replacement , ")" , collapse="" )
        
        w.start <- substr(s,1,1)
        w.end <- substr(s,nchar(s),nchar(s))
        
        r <- paste( w.start , replacement , w.end , sep="" , collapse="" )
        gsub( pattern=s , replacement=r , x=x , fixed=TRUE )
    }
    
    ########################################
    # convert formulas to text in proper call order
    
    # check for vectorized priors and expand
    # check for link functions on left and convert to inv-link on right
    # also build list of parameters with priors as we go
    pars_with_priors <- list()
    if ( length(flist) > 1 ) {
        flag_flatten <- FALSE
        for ( i in 2:length(flist) ) {
            if ( !(class(flist[[i]])=="formula") )
                stop( "Input not a formula." )
            LHS <- flist[[i]][[2]]
            if ( class(LHS)=="call" ) {
                fname <- as.character(LHS[[1]])
                if ( fname=="c" | fname=="[" | fname %in% link.names ) {
                    if ( fname=="c" ) {
                        # found a vector prior
                        newflist <- list()
                        num_pars <- length(LHS) - 1
                        for ( j in 1:num_pars ) {
                            newflist[[j]] <- flist[[i]]
                            newflist[[j]][[2]] <- LHS[[j+1]]
                            pars_with_priors[[ as.character(LHS[[j+1]]) ]] <- 1
                        } #j
                        flist[[i]] <- newflist
                        flag_flatten <- TRUE
                    }
                    if ( fname %in% link.names ) {
                        the.link <- as.character( LHS[[1]] )
                        # strip link from left hand side
                        flist[[i]][[2]] <- flist[[i]][[2]][[2]]
                        # move inverse link to right hand side
                        the.invlink <- invlink.names[ which(link.names==the.link) ]
                        old.RHS <- flist[[i]][[3]]
                        flist[[i]][[3]] <- as.call( list( as.name(the.invlink) , old.RHS ) ) 
                    }
                    if ( fname=="[" ) {
                        # par[.], so need to expand any index variable
                        pars_with_priors[[deparse(LHS)]] <- 1
                    }
                } else {
                    # a call other than c() detected on left hand side
                    stop( paste( "Invalid prior specification:" , deparse(flist[[i]]) ) )
                }
            } else {
                # not a call, so just record name of parameter
                pars_with_priors[[as.character(LHS)]] <- 1
            }
        } #i
        if ( flag_flatten ) flist <- unlist(flist,recursive=FALSE)
    }
    if (debug) print(flist)
    
    # convert from formulas to text
    # this also splits linear models
    flist2 <- lapply( flist , formula2text )

    if (debug) print(flist2)
    
    links <- list()
    
    # check for linear models in flist2 and do search-replace in likelihood
    if ( length(flist2) > 1 ) {
        # loop in reverse order, so linear models lower down get pulled up
        for ( i in length(flist2):2 ) {
            # linear models are class list
            if ( class(flist2[[i]])=="list" ) {
                LHS <- flist2[[i]][[1]]
                RHS <- flist2[[i]][[2]]
                # save current likelihood, so can check for link
                lik_save <- flist2[[1]]
                # replace in likelihood
                flist2[[1]] <- mygrep( LHS , RHS , flist2[[1]] , add.par=FALSE )
                
                # need a link function?
                if ( flist2[[1]] != lik_save ) {
                    # build a link function with current linear model in it
                    links[[ length(links)+1 ]] <- flist2[[i]]
                }
                
                # also search in other linear models above this one
                if ( i > 2 ) {
                    # RHSp <- paste( "(" , RHS , ")" , collapse="" )
                    for ( j in (i-1):2 ) {
                        if ( class(flist2[[j]])=="list" ) {
                            #flist2[[j]][[2]] <- gsub( LHS , RHSp , flist2[[j]][[2]] )
                            flist2[[j]][[2]] <- mygrep( LHS , RHS , flist2[[j]][[2]] , add.par=TRUE )
                        }
                    }
                }
                # remove symbol from list of parameters with priors
                pars_with_priors[[ LHS ]] <- NULL
            }
        }
        # second pass to remove linear models
        flist3 <- list()
        j <- 1
        for ( i in 1:length(flist2) ) {
            if ( class(flist2[[i]]) != "list" ) {
                flist3[[j]] <- flist2[[i]]
                j <- j + 1
            }
        }
        flist2 <- flist3
    }

    flist.ll <- flist2[[1]] # first formula, just for log-likelihood

    if (debug) {
        print(flist.ll)
    }
    
    ########################################
    # prep and check start list
    
    if (debug) {
        print(start)
        print(pars_with_priors)
    }
    pars <- start

    # make sure all parameters with priors are in start list
    # if priors defined, but not in start list, use prior to guess start value
    if ( any( !(names(pars_with_priors) %in% names(pars)) ) ) {
        bad_pars <- names(pars_with_priors)[ !(names(pars_with_priors) %in% names(pars)) ]
        
        #stop( paste( "Priors defined for parameters not in start list:" , paste(bad_pars,collapse=" ") ) )
        
        if ( verbose==TRUE )
            message( paste( "Sampling start values from priors for:" , paste(bad_pars,collapse=" ") ) )
        
        for ( k in bad_pars ) {
            # scan formula for right prior
            for ( g in 2:length(flist) ) {
                # check for `[`
                if ( class(flist[[g]][[2]])=="call" ) {
                    # assume `[`, because other calls should be purged by now
                    the_par_with_index <- deparse(flist[[g]][[2]])
                    the_par <- as.character( flist[[g]][[2]][[2]] )
                    the_index_or_var <- as.character( flist[[g]][[2]][[3]] )
                    if ( the_par_with_index == k ) {
                        # right prior, but need to decide if numeric index or variable index
                        idx_num <- suppressWarnings( as.numeric(the_index_or_var) )
                        if ( !is.na(idx_num) ) {
                            # for par[num], 
                            # sample one value and insert into vector in start list
                            if ( is.null( start[[the_par]] ) ) {
                                # not in start list yet, so add
                                # get length of vector by scanning priors
                                max_index <- 1
                                for ( h in 2:length(flist) ) {
                                    if ( class(flist[[h]][[2]])=="call" ) {
                                        if ( as.character(flist[[h]][[2]][[2]])=="[" & as.character(flist[[h]][[2]][[3]])==the_par ) {
                                            nval <- suppressWarnings( as.numeric(flist[[h]][[2]][[3]]) )
                                            if ( !is.null(nval) ) {
                                                max_index <- max( max_index , nval )
                                            }#is numeric
                                        }#is `[` and has same parameter name
                                    }#is call
                                }#h
                                start[[the_par]] <- rep(0,max_index)
                                the_index <- as.numeric( flist[[g]][[2]][[3]] )
                                start[[the_par]][the_index] <- sample_from_prior( flist[[g]] )
                            } else {
                                # vector already there, so insert in right place
                                the_index <- as.numeric( flist[[g]][[2]][[3]] )
                                start[[the_par]][the_index] <- sample_from_prior( flist[[g]] )
                            }
                        }# numeric index
                        if ( is.na(idx_num) ) {
                            # par[var] format
                            # need to sample vector of values
                            # get length of index variable
                            the_var <- as.character( flist[[g]][[2]][[3]] )
                            n_unique <- length(unique(data[[the_var]]))
                            start[[the_par]] <- replicate( n_unique , sample_from_prior( flist[[g]] ) )
                        }
                    }#matched
                    
                } else {
                    # ordinary parameter name, we hope
                    the_par <- paste( as.character(flist[[g]][[2]]) , collapse="" )
                    if ( the_par == k ) {
                        f <- flist[[g]]
                        theta <- sample_from_prior( f )
                        start[[k]] <- theta
                    }
                }
            }#g
        }#k
    }
    
    # list parameters without explicit priors and warn
    if ( FALSE ) {
        if ( any( !(names(pars) %in% names(pars_with_priors)) ) ) {
            flat_pars <- names(pars)[ !(names(pars) %in% names(pars_with_priors)) ]
            message( paste( "Using flat priors for:" , paste(flat_pars,collapse=" ") ) )
        }
    }
    
    #################
    # handle vector parameters
    # at this point in parse, may have vector paramters in pars (and start)
    # optim can't handle vector paramters
    # so need to:
    # (1) convert vector to series of scalar par_._1, par_._2, ..., par_._n
    # (2) pass new pars that are flat scalar to optim
    # (3) inside minuslogl, parameters names par_._1 etc are converted back to vector with name 'par'
    # (4) after optim converges, convert solution back to vector, for both coef and vcov
    
    pars <- start
    
    pars_flat <- list()
    veclist <- list() # holds info needed to convert back to scalars during search
    for ( i in 1:length(pars) ) {
        n <- length( pars[[i]] )
        par_name <- names(pars)[i]
        if ( n ==1 ) {
            # regular parameter, so just copy it
            pars_flat[[par_name]] <- pars[[i]]
        } else {
            # vector parameter, so explode it into individual parameters
            for ( j in 1:n ) {
                new_name <- concat( par_name , idx_marker_string , j )
                pars_flat[[new_name]] <- pars[[i]][j]
            }#j
            # and store name and dimension info, so can easily explode later
            veclist[[par_name]] <- list( name=par_name , n=n )
        }
    }#i
    pars <- pars_flat
    
    ########################################
    # call optim for search
    fit <- try(
        suppressWarnings(optim( par=pars , fn=make_minuslogl , flist=flist2 , data=data , veclist=veclist , hessian=hessian , method=method , ... ))
        , silent=TRUE
    )
    if ( class(fit)=="try-error" ) {
        # something went wrong...try to figure it out
        msg <- attr(fit,"condition")$message
        
        objnotfound <- grep( "object '.*' not found" , msg )
        if ( length(objnotfound) > 0 ) {
            # a parameter without prior or start value?
            obj_name <- regmatches( msg , gregexpr( "'.*'" , msg ) )
            out_msg <- paste( "Cannot find ",obj_name,".\nIf this is a parameter, try defining a prior for it or providing a start value.\nIf this is a variable, make sure it is in the data list." , collapse="" , sep="" )
            stop( out_msg )
        }
        
        objnotfound <- grep( "initial value in 'vmmin' is not finite" , msg )
        if ( length(objnotfound) > 0 ) {
            # missing values in data?
            out_msg <- paste( msg,"\nThe start values for the parameters were invalid. This could be caused by missing values (NA) in the data or by start values outside the parameter constraints. If there are no NA values in the data, try using explicit start values." , collapse="" , sep="" )
            stop( out_msg )
        }
        
        objnotfound <- grep( "non-finite finite-difference value" , msg )
        if ( length(objnotfound) > 0 ) {
            # bad start values?
            out_msg <- paste( msg,"\nStart values for parameters may be too far from MAP.\nTry better priors or use explicit start values.\nIf you sampled random start values, just trying again may work.\nStart values used in this attempt:\n", paste(names(start),"=",start,collapse="\n") , collapse="" , sep="" )
            stop( out_msg )
        }
        
        # not recognized
        stop(msg)
    }
    
    ########################################
    # prep results
    
    # convert vector parameters back to vector form
    #pars_vec <- pars_to_vectors( fit$pars , veclist )
    #pars_raw <- fit$pars
    
    if ( hessian ) {
        vcov <- try( solve(fit$hessian) )
        if ( class(vcov)=="try-error" ) {
            warning( "Error when computing variance-covariance matrix (Hessian). Fit may not be reliable." )
            vcov <- matrix( NA , nrow=length(pars) , ncol=length(pars) )
        }
    } else {
        vcov <- matrix( NA , nrow=length(pars) , ncol=length(pars) )
    }
    
    # compute minus log-likelihood at MAP, ignoring priors to do so
    # need this for correct deviance calculation, as deviance ignores priors
    fit$minuslogl <- make_minuslogl( fit$par , flist=flist.ll , data=data , veclist=veclist )
    
    # function to use later in computing DIC
    fmll <- function(pars) make_minuslogl( pars , flist=flist.ll , data=data , veclist=veclist )
    
    # rename any _._n parameters to [n]
    coefs <- fit$par
    for ( i in 1:length(coefs) ) {
        a_split <- strsplit( names(coefs)[i] , idx_marker_string )[[1]]
        if ( length(a_split)>1 ) {
            new_name <- concat( a_split[1] , "[" , a_split[2] , "]" )
            names(coefs)[i] <- new_name
        }
    }
    
    ########################################
    # build and return result
    m <- new( "map" , 
            call = match.call(), 
            coef = coefs, 
            vcov = vcov,
            optim = fit,
            data = as.list(data),
            start = start,
            formula = flist.orig,
            formula_parsed = flist2,
            fminuslogl = fmll,
            links = links )
    attr(m,"df") <- length(m@coef)
    attr(m,"veclist") <- veclist
    if (!missing(data)) attr(m,"nobs") = length(data[[1]])
    # check convergence and warn
    xcheckconvergence(m)
    # result
    m
}

# EXAMPLES
if ( FALSE ) {

data(cars)

flist0 <- list(
    dist ~ dnorm( mean=a+b*speed , sd=sigma )
)

flist0 <- list(
    dist ~ dnorm( mean=mu , sd=sigma ) ,
    mu ~ a+b*speed
)

flist1 <- list(
    dist ~ dnorm( mean=a+b*speed , sd=sigma ) ,
    b ~ dnorm(0,1) ,
    sigma ~ dcauchy(0,1)
)

flist1 <- list(
    dist ~ dnorm( mean=mu , sd=sigma ) ,
    mu ~ a+b*speed ,
    b ~ dnorm(0,1) ,
    sigma ~ dcauchy(0,1)
)

flist2 <- list(
    dist ~ dnorm( mean=a+b*speed , sd=sigma ) ,
    c(a,b) ~ dnorm(0,10) , 
    sigma ~ dcauchy(0,1)
)

# curve( dlaplace(x,1,0) , from=-10, to=10 )

flist3 <- list(
    dist ~ dnorm( mean=a+b*speed , sd=sigma ) ,
    b ~ dlaplace(1) , 
    sigma ~ dcauchy(0,1)
)

fit <- map( flist1 , start=list(a=40,b=0.1,sigma=20) , data=cars , debug=FALSE )

#########

library(rethinking)
data(chimpanzees)

flist1 <- list(
    pulled.left ~ dbinom( prob=logistic( a + b*prosoc.left ) , size=1 ),
    c(a,b) ~ dnorm(0,1)
)

flist2 <- list(
    pulled.left ~ dbinom( size=1 , prob=logistic( a + b*prosoc.left ) ),
    b ~ dnorm(0,1)
)

flist4 <- alist(
    pulled.left ~ dbinom( prob=p , size=1 ),
    logit(p) <- a + b*prosoc.left ,
    c(a,b) ~ dnorm(0,1)
)

fit2 <- map( flist4 , data=chimpanzees , start=list(a=0,b=0) , debug=FALSE )

########
# regularized logistic regression example
y <- c( rep(0,10) , rep(1,10) )
x <- c( rep(-1,9) , rep(1,11) )

flist0 <- list(
    y ~ dbinom( prob=logistic( a + b*x ) , size=1 )
)

flist1 <- list(
    y ~ dbinom( prob=logistic( a + b*x ) , size=1 ),
    c(a,b) ~ dnorm(0,10)
)

fit3a <- map( flist0 , data=list(y=y,x=x) , start=list(a=0,b=0) )

fit3b <- map( flist1 , data=list(y=y,x=x) , start=list(a=0,b=0) )

plot( y ~ x )
p <- sample.naive.posterior(fit3b)
xseq <- seq(-1,1,length.out=20)
pi.mu <- sapply( xseq , function(x) mean(logistic(p$a+p$b*x)) )
pi.ci <- sapply( xseq , function(x) PCI(logistic(p$a+p$b*x)) )
lines( xseq , pi.mu )
shade( pi.ci , xseq )

} #EXAMPLES