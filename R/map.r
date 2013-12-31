# MAP - maximum a posteriori estimation

map <- function( flist , start , data , method="BFGS" , hessian=TRUE , debug=FALSE , ... ) {
    
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
    if ( missing(start) ) stop( "'start' list required." )
    if ( missing(data) ) stop( "'data' required." )
    if ( !( class(data) %in% c("list","data.frame") ) ) {
        stop( "'data' must be of class list or data.frame." )
    }
    
    flist.orig <- flist
    
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
    
    # these function are the engine that does the optimization
    # trick is to build R code for objective function that is passed to parser
    # it's slow, but easier to debug

    dparser <- function( flist , e ) {
    # merge this into make_minuslogl?
        r <- sapply( flist , function(i) sum(eval(parse(text=i),envir=e)) )
        -sum(r)
    }
    
    make_minuslogl <- function( pars , flist , data ) {
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
        fname <- as.character( RHS[[1]] )
        if ( fname=="+" | fname=="*" | fname=="-" | fname=="/" | fname=="%*%" | fname %in% invlink.names ) {
            # linear model formula with no density (but maybe invlink) function
            # return a list with parameter name in [[1]] and text of RHS in [[2]]
            thetext <- list( as.character(LHS) , paste( deparse(RHS) , collapse=" " ) )
        } else {
            # likelihood or prior formula
            n_args <- length(RHS)
            args_list <- as.list(RHS)
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
    
    mygrep <- function( target , replacement , x , add.par=TRUE ) {
        wild <- "[()=*+ ]"
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
                if ( as.character(LHS[[1]])=="c" | as.character(LHS[[1]]) %in% link.names ) {
                    if ( as.character(LHS[[1]])=="c" ) {
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
                    if ( as.character(LHS[[1]]) %in% link.names ) {
                        the.link <- as.character( LHS[[1]] )
                        # strip link from left hand side
                        flist[[i]][[2]] <- flist[[i]][[2]][[2]]
                        # move inverse link to right hand side
                        the.invlink <- invlink.names[ which(link.names==the.link) ]
                        old.RHS <- flist[[i]][[3]]
                        flist[[i]][[3]] <- as.call( list( as.name(the.invlink) , old.RHS ) ) 
                    }
                } else {
                    # a call other than c() detected on left hand side
                    stop( paste("Invalid prior specification:",flist[[i]]) )
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
    
    # check for linear models in flist2 and do search-replace in likelihood
    if ( length(flist2) > 1 ) {
        for ( i in length(flist2):2 ) {
            # linear models are class list
            if ( class(flist2[[i]])=="list" ) {
                LHS <- flist2[[i]][[1]]
                RHS <- flist2[[i]][[2]]
                # replace in likelihood
                #flist2[[1]] <- gsub( LHS , RHS , flist2[[1]] )
                flist2[[1]] <- mygrep( LHS , RHS , flist2[[1]] , add.par=FALSE )
                
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
    pars <- unlist(start)

    if (debug) {
        print(pars)
        print(pars_with_priors)
    }

    # make sure all parameters with priors are in start list
    if ( any( !(names(pars_with_priors) %in% names(pars)) ) ) {
        bad_pars <- names(pars_with_priors)[ !(names(pars_with_priors) %in% names(pars)) ]
        stop( paste( "Priors defined for parameters not in start list:" , paste(bad_pars,collapse=" ") ) )
    }
    # list parameters without explicit priors and warn
    if ( FALSE ) {
        if ( any( !(names(pars) %in% names(pars_with_priors)) ) ) {
            flat_pars <- names(pars)[ !(names(pars) %in% names(pars_with_priors)) ]
            message( paste( "Using flat priors for:" , paste(flat_pars,collapse=" ") ) )
        }
    }
    
    ########################################
    # call optim for search
    fit <- optim( par=pars , fn=make_minuslogl , flist=flist2 , data=data , hessian=hessian , method=method , ... )
    
    ########################################
    # prep results
    
    if ( hessian ) {
        vcov <- solve(fit$hessian)
    } else {
        vcov <- matrix( NA , nrow=length(pars) , ncol=length(pars) )
    }
    
    # compute minus log-likelihood at MAP, ignoring priors to do so
    # need this for correct deviance calculation, as deviance ignores priors
    fit$minuslogl <- make_minuslogl( fit$par , flist=flist.ll , data=data )

    fmll <- function(pars) make_minuslogl( pars , flist=flist.ll , data=data )
    
    ########################################
    # build and return result
    m <- new( "map" , 
            call = match.call(), 
            coef = fit$par, 
            vcov = vcov,
            optim = fit,
            data = as.list(data),
            formula = flist.orig,
            fminuslogl = fmll )
    attr(m,"df") = length(m@coef)
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

data(chimpanzees)

flist1 <- list(
    pulled.left ~ dbinom( prob=logistic( a + b*prosoc.left ) , size=1 ),
    c(a,b) ~ dnorm(0,1)
)

flist2 <- list(
    pulled.left ~ dbinom( size=1 , prob=logistic( a + b*prosoc.left ) ),
    b ~ dnorm(0,1)
)

flist4 <- list(
    pulled.left ~ dbinom( prob=p , size=1 ),
    logit(p) ~ a + b*prosoc.left ,
    c(a,b) ~ dnorm(0,1)
)

fit2 <- map( flist3 , data=chimpanzees , start=list(a=0,b=0) , debug=TRUE )

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