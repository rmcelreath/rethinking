# glimmer
# translate glmer style model formulas into equivalent formula alist suitable for map2stan
# just a stub for now

glimmer <- function( formula , family="Gaussian" , ... ) {

    # templates
    family_liks <- list(
        Gaussian = "dnorm( mu , sigma )",
        Binomial = "dbinom( size , p )",
        Poisson = "dpois( lambda )"
    )
    
    # check input
    if ( class(formula)!="formula" ) stop( "Input must be a glmer-style formula." )
    
    f <- formula
    flist <- alist()
    
    # get outcome
    outcome <- f[[2]]
    
    # build likelihood
    flist[[1]] <- as.formula( concat( as.character(outcome) , " ~ " , family_liks[[family]] ) )
    
    # build linear model
    
    
}


parse_glimmer_formula <- function( formula , data ) {
    ## take a formula and parse into fixed effect and varying effect lists
    nobars <- function(term) {
        if (!('|' %in% all.names(term))) return(term)
        if (is.call(term) && term[[1]] == as.name('|')) return(NULL)
        if (length(term) == 2) {
        nb <- nobars(term[[2]])
        if (is.null(nb)) return(NULL)
        term[[2]] <- nb
        return(term)
        }
        nb2 <- nobars(term[[2]])
        nb3 <- nobars(term[[3]])
        if (is.null(nb2)) return(nb3)
        if (is.null(nb3)) return(nb2)
        term[[2]] <- nb2
        term[[3]] <- nb3
        term
    }
    findbars <- function(term) {
        if (is.name(term) || !is.language(term)) return(NULL)
        if (term[[1]] == as.name("(")) return(findbars(term[[2]]))
        if (!is.call(term)) stop("term must be of class call")
        if (term[[1]] == as.name('|')) return(term)
        if (length(term) == 2) return(findbars(term[[2]]))
        c(findbars(term[[2]]), findbars(term[[3]]))
    }
    subbars <- function(term)
    ### Substitute the '+' function for the '|' function
    {
        if (is.name(term) || !is.language(term)) return(term)
        if (length(term) == 2) {
        term[[2]] <- subbars(term[[2]])
        return(term)
        }
        stopifnot(length(term) >= 3)
        if (is.call(term) && term[[1]] == as.name('|'))
        term[[1]] <- as.name('+')
        for (j in 2:length(term)) term[[j]] <- subbars(term[[j]])
        term
    }
    hasintercept <- function(term) {
        attr( terms(term) , "intercept" )==1
    }
    
    # find fixed effects list by deleting random effects and expanding
    f_nobars <- nobars( formula )
    # catch implied intercept error -- happens when right side of formula is only () blocks
    if ( class(f_nobars)=="name" & length(f_nobars)==1 ) {
        f_nobars <- nobars( as.formula( paste( deparse(formula) , "+ 1" ) ) )
    }
    fixef <- colnames( model.matrix( f_nobars , data ) )
    
    # convert to all fixed effects and build needed model matrix
    mdat <- model.matrix( subbars( formula ) , data )
    outcome_name <- deparse( f_nobars[[2]] )
    # mdat <- cbind( data[[outcome_name]] , mdat )
    # colnames(mdat)[1] <- outcome_name
    outcome <- model.frame( f_nobars , data )[,1]
    if ( class(outcome)=="matrix" ) {
        # fix outcome name
        outcome_name <- colnames( outcome )[1]
    }
    
    # check for any varying effects
    if ( formula == nobars(formula) ) {
        # no varying effects
        ranef <- list()
    } else {
        # find varying effects list
        var <- findbars( formula )
        ranef <- list()
        for ( i in 1:length(var) ) {
            name <- all.vars( var[[i]][[3]] )

            if ( FALSE ) { # check index variables?
                if ( TRUE ) {
                    # check that grouping vars are class integer
                    if ( class( data[[name]] )!="integer" ) {
                        stop( paste( "Grouping variables must be integer type. '" , name , "' is instead of type: " , class( data[[name]] ) , "." , sep="" ) )
                    }
                    # check that values are contiguous
                    if ( min(data[[name]]) != 1 ) 
                        stop( paste( "Group variable '" , name , "' doesn't start at index 1." , sep="" ) )
                    ulist <- unique( data[[name]] )
                    diffs <- ulist[2:length(ulist)] - ulist[1:(length(ulist)-1)]
                    if ( any(diffs != 1) )
                        stop( paste( "Group variable '" , name , "' is not contiguous." , sep="" ) )
                } else {
                    # new code to coerce grouping vars to type integer
                    # first, save the labels by converting to character
                    #glabels <- as.character( data[[name]] )
                    # second, coerce to factor then integer to make contiguous integer index
                    mdat[,name] <- as.integer( as.factor( mdat[,name] ) )
                }
            }
            
            # parse formula
            v <- var[[i]][[2]]
            if ( class(v)=="numeric" ) {
                # just intercept
                ranef[[ name ]] <- "(Intercept)"
            } else {
                # should be class "call" or "name"
                # need to convert formula to model matrix headers
                f <- as.formula( paste( "~" , deparse( v ) , sep="" ) )
                ranef[[ name ]] <- colnames( model.matrix( f , data ) )
            }
        }
    }
    
    # result sufficient information to build Stan model code
    list( y=outcome , yname=outcome_name , fixef=fixef , ranef=ranef , dat=as.data.frame(mdat) )
}

