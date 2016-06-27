# glimmer
# translate glmer style model formulas into equivalent formula alist suitable for map2stan
# just a stub for now

xparse_glimmer_formula <- function( formula , data ) {
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
    #fixef <- make.names( colnames( model.matrix( f_nobars , data ) ) )
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


glimmer <- function( formula , data , family=gaussian , prefix=c("b_","v_") , default_prior="dnorm(0,10)" , ... ) {
    
    undot <- function( astring ) {
        astring <- gsub( "." , "_" , astring , fixed=TRUE )
        astring <- gsub( ":" , "_X_" , astring , fixed=TRUE )
        astring <- gsub( "(" , "" , astring , fixed=TRUE )
        astring <- gsub( ")" , "" , astring , fixed=TRUE )
        astring
    }
    
    # convert family to text
    family.orig <- family
    if ( class(family)=="function" ) {
        family <- do.call(family,args=list())
    }
    link <- family$link
    family <- family$family
    
    # templates
    family_liks <- list(
        gaussian = "dnorm( mu , sigma )",
        binomial = "dbinom( size , p )",
        poisson = "dpois( lambda )"
    )
    lm_names <- list(
        gaussian = "mu",
        binomial = "p",
        poisson = "lambda"
    )
    link_names <- list(
        gaussian = "identity",
        binomial = "logit",
        poisson = "log"
    )
    
    # check input
    if ( class(formula)!="formula" ) stop( "Input must be a glmer-style formula." )
    if ( missing(data) ) stop( "Need data" )
    
    f <- formula
    flist <- alist()
    prior_list <- alist()
    
    # parse
    pf <- xparse_glimmer_formula( formula , data )
    pf$yname <- undot(pf$yname)
    
    # build likelihood
    # check for size variable in Binomial
    dtext <- family_liks[[family]]
    if ( family=="binomial" ) {
        if ( class(pf$y)=="matrix" ) {
            # cbind input
            pf$dat[[pf$yname]] <- pf$y[,1]
            pf$dat[[concat(pf$yname,"_size")]] <- apply(pf$y,1,sum)
            dtext <- concat("dbinom( ",concat(pf$yname,"_size")," , p )")
        } else {
            # bernoulli
            pf$dat[[pf$yname]] <- pf$y
            dtext <- concat("dbinom( 1 , p )")
        }
    } else {
        pf$dat[[pf$yname]] <- pf$y
    }
    flist[[1]] <- concat( as.character(pf$yname) , " ~ " , dtext )
    
    # build fixed linear model
    flm <- ""
    for ( i in 1:length(pf$fixef) ) {
        # for each term, add corresponding term to linear model text
        aterm <- undot(pf$fixef[i])
        newterm <- ""
        if ( aterm=="Intercept" ) {
            newterm <- aterm
            prior_list[[newterm]] <- default_prior
        } else {
            par_name <- concat( prefix[1] , aterm )
            newterm <- concat( par_name , "*" , aterm )
            prior_list[[par_name]] <- default_prior
        }
        if ( i > 1 ) flm <- concat( flm , " +\n        " )
        flm <- concat( flm , newterm )
    }
    
    vlm <- ""
    num_group_vars <- length(pf$ranef)
    if ( num_group_vars > 0 ) {
    for ( i in 1:num_group_vars ) {
        group_var <- undot(names(pf$ranef)[i])
        members <- list()
        for ( j in 1:length(pf$ranef[[i]]) ) {
            aterm <- undot(pf$ranef[[i]][j])
            newterm <- ""
            var_prefix <- prefix[2]
            if ( num_group_vars>1 ) var_prefix <- concat( var_prefix , group_var , "_" )
            if ( aterm=="Intercept" ) {
                par_name <- concat( var_prefix , aterm )
                newterm <- concat( par_name , "[" , group_var , "]" )
            } else {
                par_name <- concat( var_prefix , aterm )
                newterm <- concat( par_name , "[" , group_var , "]" , "*" , aterm )
            }
            members[[par_name]] <- par_name
            if ( i > 1 | j > 1 ) vlm <- concat( vlm , " +\n        " )
            vlm <- concat( vlm , newterm )
        }#j
        # add group prior
        if ( length(members)>1 ) {
            # multi_normal
            gvar_name <- concat( "c(" , paste(members,collapse=",") , ")" , "[" , group_var , "]" )
            prior_list[[gvar_name]] <- concat( "dmvnorm2(0,sigma_" , group_var , ",Rho_" , group_var , ")" )
            prior_list[[concat("sigma_",group_var)]] <- concat( "dcauchy(0,2)" )
            prior_list[[concat("Rho_",group_var)]] <- concat( "dlkjcorr(2)" )
        } else {
            # normal
            gvar_name <- concat( members[[1]] , "[" , group_var , "]" )
            prior_list[[gvar_name]] <- concat( "dnorm(0,sigma_" , group_var , ")" )
            prior_list[[concat("sigma_",group_var)]] <- concat( "dcauchy(0,2)" )
        }
        # make sure grouping variables in dat
        # also ensure is an integer index
        pf$dat[[group_var]] <- coerce_index( data[[group_var]] )
    }#i
    }# ranef processing
    
    # any special priors for likelihood function
    if ( family=="gaussian" ) {
        prior_list[["sigma"]] <- "dcauchy(0,2)"
    }
    
    # insert linear model
    lm_name <- lm_names[[family]]
    #link_func <- link_names[[family]]
    if ( vlm=="" )
        lm_txt <- concat( flm )
    else
        lm_txt <- concat( flm , " +\n        " , vlm )
    lm_left <- concat( link , "(" , lm_name , ")" )
    if ( link=="identity" ) lm_left <- lm_name
    flist[[2]] <- concat( lm_left , " <- " , lm_txt )
    
    # build priors
    for ( i in 1:length(prior_list) ) {
        pname <- names(prior_list)[i]
        p_txt <- prior_list[[i]]
        flist[[i+2]] <- concat( pname , " ~ " , p_txt )
    }
    
    # build formula text and parse into alist object
    flist_txt <- "alist(\n"
    for ( i in 1:length(flist) ) {
        flist_txt <- concat( flist_txt , "    " , flist[[i]] )
        if ( i < length(flist) ) flist_txt <- concat( flist_txt , ",\n" )
    }
    flist_txt <- concat( flist_txt , "\n)" )
    flist2 <- eval(parse(text=flist_txt))
    
    # clean variable names
    names(pf$dat) <- sapply( names(pf$dat) , undot )
    # remove Intercept from dat
    pf$dat[['Intercept']] <- NULL
    
    # result
    cat(flist_txt)
    cat("\n")
    invisible(list(f=flist2,d=pf$dat))
    
}


####### TEST CODE ########
if ( FALSE ) {

library(rethinking)

data(chimpanzees)
f0 <- pulled.left ~ prosoc.left*condition - condition
m0 <- glimmer( f0 , chimpanzees , family=binomial )

f1 <- pulled.left ~ (1|actor) + prosoc.left*condition - condition
m1 <- glimmer( f1 , chimpanzees , family=binomial )
# m1s <- map2stan( m1$f , data=m1$d , sample=TRUE )

f2 <- pulled.left ~ (1+prosoc.left|actor) + prosoc.left*condition - condition
m2 <- glimmer( f2 , chimpanzees , family=binomial )

data(UCBadmit)
f3 <- cbind(admit,reject) ~ (1|dept) + applicant.gender
m3 <- glimmer( f3 , UCBadmit , binomial )
m3s <- map2stan( m3$f , data=m3$d )

f4 <- cbind(admit,reject) ~ (1+applicant.gender|dept) + applicant.gender
m4 <- glimmer( f4 , UCBadmit , binomial )
m4s <- map2stan( m4$f , data=m4$d )

data(Trolley)
f5 <- response ~ (1|id) + (1|story) + action + intention + contact
m5 <- glimmer( f5 , Trolley )
m5s <- map2stan( m5$f , m5$d , sample=FALSE )

f6 <- response ~ (1+action+intention|id) + (1+action+intention|story) + action + intention + contact
m6 <- glimmer( f6 , Trolley )
m6s <- map2stan( m6$f , m6$d , sample=TRUE )

}

