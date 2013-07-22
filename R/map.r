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
    
    formula2text <- function( f ) {
    # must keep named arguments with names, 
    # so user can pass arguments out of order to density function
        RHS <- f[[3]]
        LHS <- f[[2]]
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
        paste( RHS[[1]] , "(" , args_text , ")" , sep="" )
    }
    
    ########################################
    # convert formulas to text in proper call order
    
    # check for vectorized priors and expand
    # also build list of parameters with priors as we go
    pars_with_priors <- list()
    if ( length(flist) > 1 ) {
        flag_flatten <- FALSE
        for ( i in 2:length(flist) ) {
            if ( !(class(flist[[i]])=="formula") )
                stop( "Input not a formula." )
            LHS <- flist[[i]][[2]]
            if ( class(LHS)=="call" ) {
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
    flist2 <- lapply( flist , formula2text )
    if (debug) print(flist2)
    
    ########################################
    # prep and check start list
    pars <- unlist(start)
    if (debug) print(pars_with_priors)
    # make sure all parameters with priors are in start list
    if ( any( !(names(pars_with_priors) %in% names(pars)) ) ) {
        bad_pars <- names(pars_with_priors)[ !(names(pars_with_priors) %in% names(pars)) ]
        stop( paste( "Priors defined for parameters not in start list:" , paste(bad_pars,collapse=" ") ) )
    }
    # list parameters without explicit priors and warn
    if ( any( !(names(pars) %in% names(pars_with_priors)) ) ) {
        flat_pars <- names(pars)[ !(names(pars) %in% names(pars_with_priors)) ]
        message( paste( "Using flat priors for:" , paste(flat_pars,collapse=" ") ) )
    }
    
    ########################################
    # call optim for search
    fit <- optim( par=pars , fn=make_minuslogl , flist=flist2 , data=data , hessian=hessian , method=method , ... )
    
    ########################################
    # prep results
    #fit$minuslogl <- make_minuslogl
    if ( hessian ) {
        vcov <- solve(fit$hessian)
    } else {
        vcov <- matrix( NA , nrow=length(pars) , ncol=length(pars) )
    }
    
    ########################################
    # build and return result
    m <- new("map" , 
            call = match.call(), 
            coef = fit$par, 
            vcov = vcov,
            optim = fit,
            data = as.list(data),
            formula = flist.orig)
    attr(m,"df") = length(m@coef)
    if (!missing(data)) attr(m,"nobs") = length(data[[1]])
    m
}

# EXAMPLES
if ( FALSE ) {

data(cars)

flist0 <- list(
    dist ~ dnorm( mean=a+b*speed , sd=sigma )
)

flist1 <- list(
    dist ~ dnorm( mean=a+b*speed , sd=sigma ) ,
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

fit <- map( flist3 , start=list(a=40,b=0.1,sigma=20) , data=cars )

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

fit2 <- map( flist2 , data=chimpanzees , start=list(a=0,b=0) , debug=FALSE )

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

p <- sample.naive.posterior(fit3b)
dens(p$b,adj=1)
curve( dnorm(x,0,10) , from=-10, to=15 , lty=2 , add=TRUE )

} #EXAMPLES