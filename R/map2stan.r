# map2stan2 - rewrite of compilation algorithm
# use templates this time
# build all parts of stan code at same time, as pass through formulas

# template design:
# templates map R density functions onto Stan functions

# to-do:
# (*) need to do p*theta and (1-p)*theta multiplication outside likelihood in betabinomial (similar story for gammapoisson) --- or is there an operator for pairwise vector multiplication?
# (-) handle improper input more gracefully
# (-) add "as is" formula type, with quoted text on RHS to dump into Stan code

##################
# map2stan itself

map2stan <- function( flist , data , start , pars , constraints=list() , types=list() , sample=TRUE , iter=2000 , warmup=floor(iter/2) , chains=1 , debug=FALSE , verbose=FALSE , WAIC=TRUE , cores=1 , rng_seed , ... ) {
    
    ########################################
    # empty Stan code
    m_data <- "" # data
    m_pars <- "" # parameters
    m_tpars1 <- "" # transformed parameters, declarations
    m_tpars2 <- "" # transformed parameters, transformation code
    m_model_declare <- "" # local declarations for linear models
    m_model_priors <- "" # all parameter sampling, incl. varying effects
    m_model_lm <- "" # linear model loops
    m_model_lik <- "" # likelihood statements at bottom
    m_model_txt <- "" # general order-trusting code
    m_gq <- "" # generated quantities, can build mostly from m_model pieces
    
    ########################################
    # inverse link list
    inverse_links <- list(
        log = 'exp',
        logit = 'inv_logit',
        logodds = 'inv_logit'
    )
    
    suffix_merge <- "_merge"
    
    templates <- map2stan.templates
    
    ########################################
    # check arguments
    if ( missing(flist) ) stop( "Formula or previous fit required." )

    if ( class(flist) == "map" ) {
        # previous map fit, so pull out components
        if ( verbose==TRUE ) message( "Using formula from map fit" )
        ftemp <- flist@formula
        if ( missing(data) ) data <- flist@data
        flist <- ftemp
    }

    previous_stanfit <- NULL
    if ( class(flist) == "map2stan" ) {
        # previous map2stan fit, so pull out components
        if ( verbose==TRUE ) message( "Using formula from map2stan fit" )
        ftemp <- flist@formula
        if ( missing(data) ) data <- flist@data
        previous_stanfit <- flist@stanfit
        flist <- ftemp
    }

    if ( class(flist) != "list" ) {
        if ( class(flist)=="formula" ) {
            flist <- list(flist)
        } else {
            if ( class(flist)!="map2stan" )
                stop( "Formula or previous map2stan fit required." )
        }
    }
    if ( missing(data) ) stop( "'data' required." )
    if ( !( class(data) %in% c("list","data.frame") ) ) {
        stop( "'data' must be of class list or data.frame." )
    }
    
    flist.orig <- flist
    flist <- flist_untag(flist) # converts <- to ~ and evals to formulas
    
    if ( missing(start) ) start <- list()
    start.orig <- start
    transpars <- list() # holds transform parameter code
    pars_elect <- list() # holds transformed pars to return in samples
    pars_hide <- list() # holds pars to hide in samples
    
    if ( iter <= warmup ) {
        # user probably meant iter as number of samples
        # so warn and convert
        warning(concat("'iter' less than or equal to 'warmup'. Setting 'iter' to sum of 'iter' and 'warmup' instead (",iter+warmup,")."))
        iter <- iter + warmup
    }
    
    # check data for scale attributes and remove
    for ( i in 1:length(data) ) {
        if ( !is.null( attr(data[[i]],"scaled:center") ) | !is.null( attr(data[[i]],"scaled:scale") ) ) {
            warning( concat("Stripping scale attributes from variable ",names(data)[i] ) )
            data[[i]] <- as.numeric(data[[i]])
        }
    }
    
    ########################################
    # private functions
    
    # this is now in rethinking namespace, so should prob remove
    concat <- function( ... ) {
        paste( ... , collapse="" , sep="" )
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
    
    wildpatt <- "[()=*+/ ,\\^]"
    
    mygrep <- function( target , replacement , x , add.par=TRUE , fixed=FALSE , wild=wildpatt ) {
        #wild <- wildpatt
        pattern <- paste( wild , target , wild , sep="" , collapse="" )
        
        m <- gregexpr( pattern , x , fixed=fixed )
        
        if ( length(m[[1]])==1 )
            if ( m==-1 ) return( x )
        
        s <- regmatches( x=x , m=m )
        
        if ( add.par==TRUE ) replacement <- paste( "(" , replacement , ")" , collapse="" )
        
        if ( class(s)=="list" ) s <- s[[1]]
        w.start <- substr(s,1,1)
        w.end <- substr(s,nchar(s),nchar(s))
        
        for ( i in 1:length(s) ) {
            r <- paste( w.start[i] , replacement , w.end[i] , sep="" , collapse="" )
            x <- gsub( pattern=s[i] , replacement=r , x=x , fixed=TRUE )
        }
        return(x)
    }
    # mygrep( "[id]" , "id[i]" , "a + b[id]" , FALSE , fixed=TRUE , wild="" )
    
    mygrep_old <- function( target , replacement , x , add.par=TRUE ) {
        wild <- wildpatt
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
    
    # for converting characters not allowed by Stan
    undot <- function( astring ) {
        astring <- gsub( "." , "_" , astring , fixed=TRUE )
        astring
    }
    
    # function for adding '[i]' to variables in linear models
    indicize <- function( target , index , x , replace=NULL ) {
        # add space to beginning and end as search buffer
        x <- paste( " " , x , " " , collapse="" , sep="" )
        y <- paste( target , "[" , index , "]" , collapse="" , sep="" )
        if ( !is.null(replace) ) {
            y <- paste( replace , "[" , index , "]" , collapse="" , sep="" )
        }
        o <- mygrep( target , y , x , FALSE )
        # remove space buffer
        substr( o , start=2 , stop=nchar(o)-1 )
    }
    
    # function for adding '[i]' to index variables already in brackets in linear models
    index_indicize <- function( target , index , x , replace=NULL ) {
        # add space to beginning and end as search buffer
        x <- paste( " " , x , " " , collapse="" , sep="" )
        target <- concat("[",target,"]")
        y <- paste( target , "[" , index , "]" , collapse="" , sep="" )
        if ( !is.null(replace) ) {
            y <- paste( replace , "[" , index , "]" , collapse="" , sep="" )
        }
        o <- mygrep( target , y , x , FALSE , fixed=TRUE , wild="" )
        # remove space buffer
        substr( o , start=2 , stop=nchar(o)-1 )
    }
    # index_indicize( "id" , "i" , "a + b[id] + g[id] + id" , replace="id" )
    
    detectvar <- function( target , x ) {
        wild <- wildpatt
        x2 <- paste( " " , x , " " , collapse="" , sep="" )
        pattern <- paste( wild , target , wild , sep="" , collapse="" )
        m <- regexpr( pattern , x2 )
        return(m)
    }
    
    # for detecting index variables already in brackets '[id]'
    detectindexvar <- function( target , x ) {
        #wild <- wildpatt
        x2 <- paste( " " , x , " " , collapse="" , sep="" )
        pattern <- paste( "[" , target , "]" , sep="" , collapse="" )
        m <- regexpr( pattern , x2 , fixed=TRUE )
        return(m)
    }
    
    # detectindexvar("id","a + b[id]")
    
    ########################################
    # parse formulas
    
    # function to sample from prior specified with density function
    sample_from_prior <- function( the_density , par_in ) {
        the_rdensity <- the_density
        substr( the_rdensity , 1 , 1 ) <- "r"
        pars <- vector(mode="list",length=length(par_in)+1)
        pars[[1]] <- 1
        for ( i in 1:(length(pars)-1) ) {
            pars[[i+1]] <- par_in[[i]]
        }
        result <- do.call( the_rdensity , args=pars )
        return(result)
    }
        
    # extract likelihood(s)
    extract_likelihood <- function( fl ) {
        #
        # test for $ function in outcome name
        if ( class(fl[[2]])=="call" ) {
            if ( as.character(fl[[2]][[1]])=="$" ) {
                outcome <- as.character( fl[[2]][[3]] )
            } else {
                # deparse to include function --- hopefully user knows what they are doing
                outcome <- deparse( fl[[2]] )
            }
        } else {
            outcome <- as.character( fl[[2]] )
        }
        template <- get_template( as.character(fl[[3]][[1]]) )
        likelihood <- template$stan_name
        
        # check for special Stan distribution with bound link function
        nat_link <- FALSE
        if ( !is.null(template$nat_link) ) nat_link <- template$nat_link
        
        likelihood_pars <- list()
        for ( i in 1:(length(fl[[3]])-1) ) likelihood_pars[[i]] <- fl[[3]][[i+1]]
        N_cases <- length( d[[ outcome ]] )
        N_name <- "N"
        nliks <- length(fp[['likelihood']])
        if ( nliks>0 ) {
            N_name <- concat( N_name , "_" , outcome )
        }
        # result
        list( 
            outcome = undot(outcome) ,
            likelihood = likelihood ,
            template = template$name ,
            pars = likelihood_pars ,
            N_cases = N_cases ,
            N_name = N_name ,
            out_type = template$out_type ,
            nat_link = nat_link
        )
    }
    
    # extract linear model(s)
    extract_linearmodel <- function( fl ) {
        # check for link function
        use_link <- TRUE
        if ( length(fl[[2]])==1 ) {
            parameter <- as.character( fl[[2]] )
            link <- "identity"
        } else {
            # link!
            parameter <- as.character( fl[[2]][[2]] )
            link <- as.character( fl[[2]][[1]] )
            use_link <- TRUE # later detect special Stan dist with built-in link
        }
        RHS <- paste( deparse(fl[[3]]) , collapse=" " )
        # find likelihood that contains this lm
        N_name <- "N" # default
        if ( length(fp[['likelihood']]) > 0 ) {
            for ( i in 1:length(fp[['likelihood']]) ) {
                pars <- fp[['likelihood']][[i]][['pars']]
                # find par in likelihood that matches lhs of this lm
                for ( j in 1:length(pars) ) {
                    if ( parameter == pars[[j]] ) {
                        # match!
                        N_name <- fp[['likelihood']][[i]][['N_name']]
                        # check if link for this linear model matches 
                        #   nat_link for likelihood
                        if ( link != "identity" ) {
                            nat_link <- fp[['likelihood']][[i]][['nat_link']]
                            if ( nat_link != FALSE ) {
                                if ( nat_link == link ) {
                                    # use special Stan distribution
                                    # so replace likelihood name in likelihood entry
                                    lik <- fp[['likelihood']][[i]][['likelihood']]
                                    lik <- concat( lik , "_" , link )
                                    fp_new <- fp
                                    fp_new[['likelihood']][[i]][['likelihood']] <- lik
                                    assign( "fp" , fp_new , inherits=TRUE )
                                    # and erase explicit link
                                    use_link <- FALSE
                                }
                            }
                        }# not identity
                    }#match
                }#j
            }#i
        }
        # result
        list(
            parameter = parameter ,
            RHS = RHS ,
            link = link ,
            N_name = N_name ,
            use_link = use_link
        )
    }
    
    # extract varying effect prior (2nd level likelihood)
    extract_vprior <- function( fl ) {
        group_var <- as.character( fl[[2]][[3]] )
        pars_raw <- fl[[2]][[2]]
        pars <- list()
        if ( length(pars_raw) > 1 ) {
            for ( i in 1:(length(pars_raw)-1) ) {
                pars[[i]] <- as.character( pars_raw[[i+1]] )
            }
        } else {
            pars[[1]] <- as.character( pars_raw )
        }
        
        tmplt <- get_template( as.character(fl[[3]][[1]]) )
        density <- tmplt$stan_name
        # input parameters
        pars_in <- list()
        for ( i in 2:length(fl[[3]]) ) pars_in[[i-1]] <- fl[[3]][[i]]
        list(
            pars_out = pars ,
            density = density ,
            pars_in = pars_in ,
            group = undot(group_var) ,
            template = tmplt$name
        )
    }
    
    # extract ordinary prior
    extract_prior <- function( fl ) {
        # apar <- as.character( fl[[2]] )
        if ( class(fl[[2]])=="call" ) {
            apar <- fl[[2]] # preserve c(a,b) call structure
        } else {
            apar <- deparse( fl[[2]] ) # deparse handles 'par[i]' correctly
        }
        template <- get_template(as.character(fl[[3]][[1]]))
        adensity <- template$stan_name
        inpars <- list()
        n <- length( fl[[3]] ) - 1
        for ( i in 1:n ) {
            inpars[[i]] <- fl[[3]][[i+1]]
        }
        list(
            par_out = apar ,
            density = adensity ,
            pars_in = inpars ,
            group = NA ,
            template = template$name
        )
    }
    
    # function to append entry to a list
    listappend <- function(x,entry) {
        n <- length(x)
        x[[n+1]] <- entry
        x
    }
    
    # function to scan distribution templates for a match in R_name or stan_name
    # returns NA if no match, name of entry otherwise
    # when there are multiple matches (for stan_name e.g.), returns first one
    scan_templates <- function( fname ) {
        the_match <- NA
        for ( i in 1:length(templates) ) {
            R_name <- templates[[i]][['R_name']]
            stan_name <- templates[[i]][['stan_name']]
            if ( fname %in% c(R_name,stan_name) ) {
                the_match <- names(templates)[i]
                return(the_match)
            }
        }
        return(the_match)
    }
    
    get_template <- function( fname ) {
        tmpname <- scan_templates(fname)
        if ( is.na(tmpname) ) stop(concat("Distribution ",fname," not recognized."))
        return(templates[[ tmpname ]])
    }
    
    # build parsed list
    fp <- list( 
        likelihood = list() ,
        lm = list() ,
        vprior = list(),
        prior = list(),
        used_predictors = list()
    )
    fp_order <- list()
    impute_bank <- list()
    d <- as.list(data)
    
    # flag variable names that contain dots
    for ( i in 1:length(d) ) {
        raw_name <- names(d)[i]
        new_name <- undot(raw_name)
        if ( raw_name != new_name )
            message( concat("Warning: Variable '",raw_name,"' contains dots '.'.\nWill attempt to remove dots internally.") )
        # names(d)[i] <- new_name
    }
    
    ####
    # pass over formulas and extract into correct slots
    # we go from top to bottom, so can see where linear models plug in, when we find them
    
    # function to check for variable name in data,
    # and if found add to used_predictors list
    tag_var <- function( var , N_name="N" ) {
        var <- undot(concat(var))
        result <- NULL
        if ( var %in% names(d) ) {
            type <- "real"
            var_class <- class(d[[var]])
            if ( var_class=="integer" ) type <- "int"
            result <- list( var=var , N=N_name , type=type )
        }
        return(result)
    }
    
    for ( i in 1:length(flist) ) {
    
        # first, scan for truncation operator T[,] on righthand side
        # marked for now by `&` function
        T_text <- ""
        RHS <- flist[[i]][[3]]
        if ( length(RHS)>1 ) {
            if ( class(RHS[[1]])=="name" ) {
                if( as.character(RHS[[1]])=="&" ) {
                    # test for T[,]
                    if ( class(RHS[[3]])=="call" ) {
                        if( length(RHS[[3]])==4 ) {
                            if( as.character(RHS[[3]][[2]])=="T" ) {
                                # got one, so deparse to text and remove from formula before parsing rest
                                T_text <- deparse( RHS[[3]] )
                                flist[[i]][[3]] <- flist[[i]][[3]][[2]] # just density part before `&`
                            }
                        }
                    }
                }
            }
        }
        
        # test for likelihood
        # if has distribution function on RHS and LHS is variable, then likelihood
        # could also be distribution assumption for predictor variable, for imputation
        is_likelihood <- FALSE
        RHS <- flist[[i]][[3]]
        if ( length(RHS) > 1 ) {
            function_name <- as.character( RHS[[1]] )
            ftemplate <- scan_templates( function_name )
            if ( !is.na(ftemplate) ) {
                # a distribution, so check for outcome variable
                LHS <- flist[[i]][[2]]
                # have to check for function call
                if ( class(LHS)=="call" ) {
                    if ( as.character(LHS[[1]])=="$" ) {
                        # strip of data frame/list name
                        LHS <- LHS[[3]]
                    } else {
                        # another function...assume variable in second position
                        LHS <- LHS[[2]]
                    }
                }
                if ( length(LHS)==1 ) {
                    # check if symbol is a variable
                    if ( (as.character(LHS)) %in% names(d) ) {
                        is_likelihood <- TRUE
                    }
                }
            }
        }
        if ( is_likelihood==TRUE ) {
            lik <- extract_likelihood( flist[[i]] )
            lik$T_text <- T_text
            n <- length( fp[['likelihood']] )
            fp[['likelihood']][[n+1]] <- lik
            fp_order <- listappend( fp_order , list(type="likelihood",i=n+1) )
            
            # add outcome to used variables
            fp[['used_predictors']][[undot(lik$outcome)]] <- list( var=undot(lik$outcome) , N=lik$N_name , type=lik$out_type )
            
            # check for missing values
            # if found, indicates predictor variable for imputation
            if ( any(is.na(d[[undot(lik$outcome)]])) ) {
                # overwrite with top 'N' name
                fp[['used_predictors']][[undot(lik$outcome)]] <- list( var=undot(lik$outcome) , N='N' , type=lik$out_type )
                
                # build info needed to perform imputation
                var_missingness <- which(is.na(d[[undot(lik$outcome)]]))
                var_temp <- ifelse( is.na(d[[undot(lik$outcome)]]) , -999 , d[[undot(lik$outcome)]] )
                # add to impute bank
                impute_bank[[ undot(lik$outcome) ]] <- list(
                    N_miss = length(var_missingness),
                    missingness = var_missingness,
                    init = mean( d[[lik$outcome]] , na.rm=TRUE )
                )
                # add missingness to data list
                missingness_name <- concat(undot(lik$outcome),"_missingness")
                N_missing_name <- concat(lik$N_name,"_missing")
                d[[ missingness_name ]] <- var_missingness
                # flag as used
                fp[['used_predictors']][[missingness_name]] <- list( var=missingness_name , N=N_missing_name , type="int" )
                # trap for only 1 missing value and vectorization issues
                if ( length(var_missingness)==1 )
                    fp[['used_predictors']][[missingness_name]] <- list( var=missingness_name , type="index" )
                # replace variable with missing values
                d[[undot(lik$outcome)]] <- var_temp
                
                # message to user
                message(concat("Imputing ",length(var_missingness)," missing values (NA) in variable '",undot(lik$outcome),"'."))
            }
            
            # check for binomial size variable and mark used
            if ( lik$likelihood=='binomial' | lik$likelihood=='beta_binomial' ) {
                sizename <- as.character(lik$pars[[1]])
                if ( !is.null( d[[sizename]] ) ) {
                    fp[['used_predictors']][[undot(sizename)]] <- list( var=undot(sizename) , N=lik$N_name , type="int" )
                }
            }
            if ( lik$template=="ZeroInflatedBinomial" ) {
                # second paramter of zibinom is binomial size variable
                sizename <- as.character(lik$pars[[2]])
                if ( !is.null( d[[sizename]] ) ) {
                    fp[['used_predictors']][[undot(sizename)]] <- list( var=undot(sizename) , N=lik$N_name , type="int" )
                }
            }
            next
        }
        
        # test for linear model
        is_linearmodel <- FALSE
        if ( length(flist[[i]][[3]]) > 1 ) {
            fname <- as.character( flist[[i]][[3]][[1]] )
            if ( fname=="[" | fname=="+" | fname=="(" | fname=="*" | fname=="-" | fname=="/" | fname=="%*%" | fname %in% c('exp','inv_logit') ) {
                is_linearmodel <- TRUE
            }
        } else {
            # RHS is length 1, so can't be a prior or density statement
            # must be a linear model with a single symbol in it, like "p ~ a"
            is_linearmodel <- TRUE
        }
        if ( is_linearmodel==TRUE ) {
            n <- length( fp[['lm']] )
            xlm <- extract_linearmodel( flist[[i]] )
            xlm$T_text <- T_text
            fp[['lm']][[n+1]] <- xlm
            fp_order <- listappend( fp_order , list(type="lm",i=n+1) )
            next
        }
        
        # test for varying effects prior
        # can use `|` or `[` to specify varying effects (vector of parameters)
        if ( length( flist[[i]][[2]] ) == 3 ) {
            flag_vprior <- FALSE
            fname <- as.character( flist[[i]][[2]][[1]] )
            if ( fname=="|" ) flag_vprior <- TRUE
            if ( fname=="[" ) {
                # test for hard-coded index
                if ( class(flist[[i]][[2]][[3]])=="name" ) {
                    # is a symbol/name, so not hard-coded index
                    flag_vprior <- TRUE
                }
            }
            if ( flag_vprior==TRUE ) {
                n <- length( fp[['vprior']] )
                xvp <- extract_vprior( flist[[i]] )
                xvp$T_text <- T_text
                fp[['vprior']][[n+1]] <- xvp
                fp_order <- listappend( fp_order , list(type="vprior",i=n+1) )
                # make sure index variabe is class 'integer'
                # if not, will be coerced and index automatically generated
                if ( !is.null(d[[ xvp$group ]]) )
                    if ( class(d[[ xvp$group ]])!="integer" )
                        d[[ xvp$group ]] <- coerce_index( d[[ xvp$group ]] )
                # next!
                next
            }
        }
        
        # an ordinary prior?
        # check for vectorized LHS
        if ( length( flist[[i]][[2]] ) > 1 ) {
            fname <- as.character( flist[[i]][[2]][[1]] )
            if ( fname=="c" ) {
                # need to distinguish between 
                    # repeat a prior, e.g.: c(a,b) ~ dnorm(0,1)
                    # multivariate prior, e.g.: c(a,b) ~ dmvnorm(0,Sigma)
                flag_mvprior <- FALSE
                if ( length(RHS)>1 ) {
                    fname <- as.character( RHS[[1]] )
                    if ( fname %in% c("dmvnorm","dmvnorm2","multi_normal") )
                        flag_mvprior <- TRUE
                }
                if ( flag_mvprior==FALSE ) {
                    # get list of parameters and add each as parsed prior
                    np <- length( flist[[i]][[2]] )
                    for ( j in 2:np ) {
                        fcopy <- flist[[i]]
                        fcopy[[2]] <- flist[[i]][[2]][[j]]
                        n <- length( fp[['prior']] )
                        xp <- extract_prior( fcopy )
                        xp$T_text <- T_text
                        fp[['prior']][[n+1]] <- xp
                        fp_order <- listappend( fp_order , list(type="prior",i=n+1) )
                    } #j
                } # not mvprior
                if ( flag_mvprior==TRUE ) {
                    # treat like ordinary prior
                    # template handles vector conversion in Stan code
                    # extract_prior() should preserve call structure
                    n <- length( fp[['prior']] )
                    xp <- extract_prior( flist[[i]] )
                    xp$T_text <- T_text
                    fp[['prior']][[n+1]] <- xp
                    fp_order <- listappend( fp_order , list(type="prior",i=n+1) )
                } # mvprior
            } # "c" call
            if ( fname=="[" ) {
                # hard-coded index
                n <- length( fp[['prior']] )
                xp <- extract_prior( flist[[i]] )
                xp$T_text <- T_text
                fp[['prior']][[n+1]] <- xp
                fp_order <- listappend( fp_order , list(type="prior",i=n+1) )
            } # "[" call
        } else {
            # ordinary simple prior
            n <- length( fp[['prior']] )
            xp <- extract_prior( flist[[i]] )
            xp$T_text <- T_text
            fp[['prior']][[n+1]] <- xp
            fp_order <- listappend( fp_order , list(type="prior",i=n+1) )
        }
    }
    
    #####
    # add index brackets in linear models
    # be careful to detect manual bracketing for parameter vectors -> only index the index variable itself
    # also need to replace names of variables 'x' with missing values for imputation with 'x_merge'
    # go through linear models
    n_lm <- length(fp[['lm']])
    if ( n_lm > 0 ) {
        for ( i in 1:n_lm ) {
            # for each variable in data, add index
            index <- "i"
            vnames <- names(d) # variables in data
            # add any linear model names, aside from this one
            if ( n_lm > 1 ) {
                for ( j in (1:n_lm)[-i] ) vnames <- c( vnames , fp[['lm']][[j]]$parameter )
            }
            for ( v in vnames ) {
                # tag if used
                used <- detectvar( v , fp[['lm']][[i]][['RHS']] )
                v_name <- v
                if ( used > -1 & v %in% names(d) ) {
                    # if variable (not lm), add to used predictors list
                    if ( is.null(fp[['used_predictors']][[undot(v)]]) ) 
                        fp[['used_predictors']][[undot(v)]] <- list( var=undot(v) , N=fp[['lm']][[i]][['N_name']] )
                    # check for imputation
                    if ( v %in% names(impute_bank) ) {
                        # use x_merge in linear model of imputed variable
                        v_name <- concat( v_name , suffix_merge )
                    }
                    if ( any(is.na(d[[v]])) ) {
                        # MISSING VALUES
                        # imputation variables should be clear by now
                        # abort compilation
                        stop( concat( "Variable '", v , "' has missing values (NA) in:\n" , fp[['lm']][[i]][['RHS']] , "\nEither remove cases with NA or declare a distribution to use for imputation." ) )
                    }
                }
                # add index and undot the name
                # do this outside the if-block above, so linear models get index too
                fp[['lm']][[i]][['RHS']] <- indicize( v , index , fp[['lm']][[i]][['RHS']] , replace=undot(v_name) )
                
                # check index variables in brackets
                used <- detectindexvar( v , fp[['lm']][[i]][['RHS']] )
                if ( used > -1 & v %in% names(d) ) {
                    # if variable (not lm), add to used predictors list
                    fp[['used_predictors']][[undot(v)]] <- list( var=undot(v) , N=fp[['lm']][[i]][['N_name']] )
                    # add index and undot the name
                    fp[['lm']][[i]][['RHS']] <- index_indicize( v , index , fp[['lm']][[i]][['RHS']] , replace=undot(v) )
                }
            }#v
            
            # for each varying effect parameter, add index with group
            # also rename parameter in linear model, so can use vector data type in Stan
            n <- length( fp[['vprior']] )
            if ( n > 0 ) {
                for ( j in 1:n ) {
                    vname <- paste( "vary_" , fp[['vprior']][[j]][['group']] , collapse="" , sep="" )
                    jindex <- paste( fp[['vprior']][[j]][['group']] , "[" , index , "]" , collapse="" , sep="" )
                    npars <- length(fp[['vprior']][[j]][['pars_out']])
                    for ( k in 1:npars ) {
                        var <- fp[['vprior']][[j]][['pars_out']][k]
                        # if only one parameter in cluster, don't need name change
                        if ( npars==1 ) {
                            fp[['lm']][[i]][['RHS']] <- indicize( var , jindex , fp[['lm']][[i]][['RHS']] )
                        } else {
                            # more than one parameter, so need vector name replacements
                            #jindexn <- paste( jindex , "," , k , collapse="" , sep="" ) 
                            jindexn <- jindex
                            #fp[['lm']][[i]][['RHS']] <- indicize( var , jindexn , fp[['lm']][[i]][['RHS']] , vname )
                            fp[['lm']][[i]][['RHS']] <- indicize( var , jindexn , fp[['lm']][[i]][['RHS']] )
                        }
                    }#k
                }#j
            }#n>0
        }#i
    }# if n_lm > 0
    
    # undot all the variable names in d
    d.orig <- d
    for ( i in 1:length(d) ) {
        oldname <- names(d)[i]
        names(d)[i] <- undot(oldname)
    }
    
    if ( debug==TRUE ) print(fp)
    
    #
    ########################################
    # build Stan code
    
    indent <- "    " # 4 spaces
    
    # holds any start values sampled from priors
    # need as many lists as chains to sample from, 
    #   so chains have different inits
    start_prior <- vector(mode="list",length=chains)
    for ( i in 1:chains ) start_prior[[i]] <- list()
    
    # pass back through parsed formulas and build Stan code
    # parsing now goes in *reverse*
    
    for ( f_num in length(fp_order):1 ) {
    
        f_current <- fp_order[[f_num]]
        
        if ( f_current$type == "prior" ) {
        
            i <- f_current$i
            
            prior <- fp[['prior']][[i]]
            tmplt <- templates[[prior$template]]
            klist <- tmplt$par_map(prior$pars_in,environment())
            for ( j in 1:length(klist) ) {
                l <- tag_var(klist[[j]])
                if ( !is.null(l) ) {
                    fp[['used_predictors']][[l$var]] <- l
                    klist[[j]] <- l$var # handles undotting
                }
            }
            if ( class(prior$par_out)!="character" )
                prior$par_out <- deparse(prior$par_out)
            parstxt <- paste( klist , collapse=" , " )
            
            # check for special formatting of left hand side
            # also handle some constraints here
            lhstxt <- prior$par_out
            if ( !is.null(tmplt$out_map) ) {
                lhstxt <- tmplt$out_map( prior$par_out , prior$pars_in , environment() )
            }
            
            # build text
            txt <- concat( indent , lhstxt , " ~ " , prior$density , "( " , parstxt , " )" , prior$T_text , ";" )
            
            #m_model_priors <- concat( m_model_priors , txt , "\n" )
            m_model_txt <- concat( m_model_txt , txt , "\n" )
            
            # check for explicit start value
            # need to clean any [] index on the parameter name
            # so use regular expression to split at '['
            par_out_clean <- strsplit( prior$par_out , "\\[" )[[1]][1]
            
            # now check if in start list
            if ( !( par_out_clean %in% names(start) ) ) {
                
                # try to get dimension of parameter from any corresponding vprior
                ndims <- 0
                if ( length(fp[['vprior']]) > 0 ) {
                    for ( vpi in 1:length(fp[['vprior']]) ) {
                        # look for parameter name in input parameters to mvnorm
                        if ( par_out_clean %in% fp[['vprior']][[vpi]]$pars_in ) {
                            ndims <- length( fp[['vprior']][[vpi]]$pars_out )
                        }
                    }
                }#find ndims
                
                # lkj_corr?
                if ( tmplt$R_name=="dlkjcorr" ) {
                    # just use identity matrix as initial value
                    if ( ndims > 0 ) {
                        for ( ch in 1:chains ) 
                            start_prior[[ch]][[ prior$par_out ]] <- diag(ndims)
                        if ( verbose==TRUE )
                            message( paste(prior$par_out,": using identity matrix as start value [",ndims,"]") )
                    } else {
                        # no vpriors parsed, so not sure what to do about lkj_corr dims
                        stop( paste(prior$par_out,": no dimension found for this matrix. Use explicit start value to define its dimension.\nFor example:",prior$par_out,"= diag(2)") )
                    }
                } else {
                    # any old anonymous prior -- sample init from prior density
                    # but must check for dimension, in case is a vector
                    ndims <- max(ndims,1)
                    theta <- list()
                    for ( ch in 1:chains ) {
                        theta[[ch]] <- try(
                            replicate( 
                                ndims , 
                                sample_from_prior( tmplt$R_name , prior$pars_in ) ),
                            silent=TRUE)
                    }#ch
                    if ( class(theta)=="try-error" ) {
                        # something went wrong
                        msg <- attr(theta,"condition")$message
                        message(concat("Error trying to sample start value for parameter '",prior$par_out,"'.\nPrior malformed, or maybe you mistyped a variable name?"))
                        stop(msg)
                    } else {
                        for ( ch in 1:chains )
                            start_prior[[ch]][[ prior$par_out ]] <- theta[[ch]]
                        if ( ndims==1 ) {
                            if ( verbose==TRUE )
                                message( paste(prior$par_out,": using prior to sample start value") )
                        } else {
                            if ( verbose==TRUE ) 
                                message( paste(prior$par_out,": using prior to sample start values [",ndims,"]") )
                        }
                    }#no error
                }
            }#not in start list
            
        } # prior
    
        if ( f_current$type == "vprior" ) {
            
            i <- f_current$i
            
            vprior <- fp[['vprior']][[i]]
            tmplt <- templates[[vprior$template]]
            N_txt <- concat( "N_" , vprior$group )
            npars <- length(vprior$pars_out)
            
            # add N count to data
            N <- length( unique( d[[ vprior$group ]] ) )
            d[[ N_txt ]] <- N
            
            # check for explicit start value
            # do this now, so out_map and pars_map functions can modify
            for ( k in vprior$pars_out )
                if ( !( k %in% names(start) ) ) {
                    if ( verbose==TRUE )
                        message( paste(k,": start values set to zero [",N,"]") )
                    for ( ch in 1:chains )
                        start_prior[[ch]][[ k ]] <- rep(0,N)
                }
            
            # lhs -- if vector, need transformed parameter of vector type
            lhstxt <- ""
            if ( length(vprior$pars_out) > 1 ) {
                # check for special formatting function for density
                if ( !is.null(tmplt$out_map) ) {
                    lhstxt <- tmplt$out_map( vprior$pars_out , vprior$pars_in , N , N_txt , environment() )
                } else {
                # automated conversion to vector
                    # parameter vector
                    lhstxt <- paste( vprior$pars_out , collapse="" )
                    lhstxt <- concat( "v_" , lhstxt )
                    
                    # add declaration to transformed parameters
                    m_tpars1 <- concat( m_tpars1 , indent , "vector[" , npars , "] " , lhstxt , "[" , N_txt , "];\n" )
                    
                    # add conversion for each true parameter
                    m_tpars2 <- concat( m_tpars2 , indent , "for ( j in 1:" , N_txt , " ) {\n" )
                    for ( j in 1:npars ) {
                        m_tpars2 <- concat( m_tpars2 , indent,indent , lhstxt , "[j," , j , "] <- " , vprior$pars_out[[j]] , "[j];\n" )
                    }
                    m_tpars2 <- concat( m_tpars2 , indent , "}\n" )
                }
            } else {
                # single parameter
                lhstxt <- vprior$pars_out[[1]]
            }
            
            # format parmater inputs
            # use par_map function in template, so ordering etc can change
            #print(length(vprior$pars_out))
            klist <- tmplt$par_map( vprior$pars_in , environment() , length(vprior$pars_out) , N_txt )
            for ( j in 1:length(klist) ) {
                l <- tag_var(klist[[j]])
                if ( !is.null(l) ) {
                    fp[['used_predictors']][[l$var]] <- l
                    klist[[j]] <- l$var # handles undotting
                }
            }
            rhstxt <- paste( klist , collapse=" , " )
            
            # add text to model code
            # new vectorized normal and multi_normal
            if ( vprior$density %in% c("multi_normal","normal") )
                m_model_txt <- concat( m_model_txt , indent , lhstxt , " ~ " , vprior$density , "( " , rhstxt , " )" , vprior$T_text , ";\n" )
            else
            # old un-vectorized code
                m_model_txt <- concat( m_model_txt , indent , "for ( j in 1:" , N_txt , " ) " , lhstxt , "[j] ~ " , vprior$density , "( " , rhstxt , " )" , vprior$T_text , ";\n" )
            
            # declare each parameter with correct type from template
            outtype <- "vector"
            for ( j in 1:length(vprior$pars_out) ) {
                #m_pars <- concat( m_pars , indent , outtype , "[" , N_txt , "] " , vprior$pars_out[[j]] , ";\n" )
            }
            
            # add data declaration for grouping variable number of unique values
            #m_data <- concat( m_data , indent , "int<lower=1> " , N_txt , ";\n" )
            fp[['used_predictors']][[N_txt]] <- list( var=N_txt , type="index" )
            
            # mark grouping variable used
            # check if already there
            if ( is.null(fp[['used_predictors']][[vprior$group]]) ) {
                # check likelihoods for matching length and use that N_name
                if ( length(fp[['likelihood']])>0 ) {
                    groupN <- length(d[[vprior$group]])
                    for ( j in 1:length(fp[['likelihood']]) ) {
                        if ( fp[['likelihood']][[j]]$N_cases == groupN ) {
                            fp[['used_predictors']][[vprior$group]] <- list( var=vprior$group , N=fp[['likelihood']][[j]]$N_name )
                            break
                        }
                    }#j
                } else {
                    # just add raw integer length
                    fp[['used_predictors']][[vprior$group]] <- list( var=vprior$group , N=length(d[[vprior$group]]) )
                }
            }# is.null
            
        } # vprior
    
        # linear models
        if ( f_current$type == "lm" ) {
        
            i <- f_current$i
        
            linmod <- fp[['lm']][[i]]
            N_txt <- linmod$N_name
            
            # if linked outcome is imputated, fall back to basic 'N'
            N_stem <- substr(N_txt,3,nchar(N_txt))
            if ( N_stem %in% names(impute_bank) ) {
                N_txt <- 'N'
            }
            
            # open the loop
            txt1 <- concat( indent , "for ( i in 1:" , N_txt , " ) {\n" )
            m_model_txt <- concat( m_model_txt , txt1 )
            m_gq <- concat( m_gq , txt1 )
            
            # assignment
            txt1 <- concat( indent,indent , linmod$parameter , "[i] <- " , linmod$RHS , ";\n" )
            m_model_txt <- concat( m_model_txt , txt1 )
            m_gq <- concat( m_gq , txt1 )
            
            # link function
            if ( linmod$link != "identity" & linmod$use_link==TRUE ) {
                # check for valid link function
                if ( is.null( inverse_links[[linmod$link]] ) ) {
                    stop( paste("Link function '",linmod$link,"' not recognized in formula line:\n",deparse(flist[[f_num]]),sep="") )
                } else {
                    # build it
                    txt1 <- concat( indent,indent , linmod$parameter , "[i] <- " , inverse_links[[linmod$link]] , "(" , linmod$parameter , "[i]);\n" )
                    m_model_txt <- concat( m_model_txt , txt1 )
                    m_gq <- concat( m_gq , txt1 )
                }
            }
            
            # close the loop
            m_model_txt <- concat( m_model_txt , indent , "}\n" )
            m_gq <- concat( m_gq , indent , "}\n" )
            
            # add declaration of linear model local variable
            # generated quantities reuse this later on in composition
            m_model_declare <- concat( m_model_declare , indent , "vector[" , N_txt , "] " , linmod$parameter , ";\n" )
            
        } # lm
    
        # likelihoods and gq
        if ( f_current$type == "likelihood" ) {
            
            i <- f_current$i

            lik <- fp[['likelihood']][[i]]
            tmplt <- templates[[lik$template]]
            parstxt_L <- tmplt$par_map( lik$pars , environment() )
            for ( j in 1:length(parstxt_L) ) {
                l <- tag_var(parstxt_L[[j]])
                if ( !is.null(l) ) {
                    fp[['used_predictors']][[l$var]] <- l
                    parstxt_L[[j]] <- l$var # handles undotting
                }
            }
            
            # add sampling statement to model block
            outcome <- lik$outcome
            
            if ( tmplt$vectorized==FALSE ) {
                # add loop for non-vectorized distribution
                txt1 <- concat( indent , "for ( i in 1:" , lik$N_name , " )\n" )
                m_model_txt <- concat( m_model_txt , txt1 )
                m_gq <- concat( m_gq , txt1 )
                
                # add [i] to outcome
                outcome <- concat( outcome , "[i]" )
                
                # check for linear model names as parameters and add [i] to each
                if ( length(fp[['lm']])>1 ) {
                    lm_names <- c()
                    for ( j in 1:length(fp[['lm']]) ) {
                        lm_names <- c( lm_names , fp[['lm']][[j]]$parameter )
                    }
                    for ( j in 1:length(parstxt_L) ) {
                        if ( as.character(parstxt_L[[j]]) %in% lm_names ) {
                            parstxt_L[[j]] <- concat( as.character(parstxt_L[[j]]) , "[i]" )
                        }
                    }
                }
            }
            
            parstxt <- paste( parstxt_L , collapse=" , " )
            
            if ( lik$likelihood=="increment_log_prob" ) {
                
                # custom distribution using increment_log_prob
                code_model <- tmplt$stan_code
                code_gq <- tmplt$stan_dev
                # replace OUTCOME, PARx symbols with actual names
                code_model <- gsub( "OUTCOME" , outcome , code_model , fixed=TRUE )
                code_gq <- gsub( "OUTCOME" , outcome , code_gq , fixed=TRUE )
                for ( j in 1:length(parstxt_L) ) {
                    parname <- as.character(parstxt_L[[j]])
                    parpat <- concat( "PAR" , j )
                    code_model <- gsub( parpat , parname , code_model , fixed=TRUE )
                    code_gq <- gsub( parpat , parname , code_gq , fixed=TRUE )
                }
                # insert into Stan code
                m_model_txt <- concat( m_model_txt , indent , code_model , "\n" )
                m_gq <- concat( m_gq , indent , code_gq , "\n" )
                
            } else {
            
                # regular distribution with ~
                
                # imputation stuff
                if ( outcome %in% names(impute_bank) ) {
                    N_miss <- impute_bank[[lik$outcome]]$N_miss
                    # add _merge suffix
                    outcome <- concat(outcome,suffix_merge)
                    # build code in transformed parameters that constructs x_merge
                    m_tpars1 <- concat( m_tpars1 , indent , "real " , outcome , "[N];\n" )
                    m_tpars2 <- concat( m_tpars2 , indent , outcome , " <- " , lik$outcome , ";\n" )
                    # trap for single missing value and vectorization issues
                    if ( N_miss>1 )
                        m_tpars2 <- concat( m_tpars2 , indent , "for ( u in 1:" , lik$N_name , "_missing ) " , outcome , "[" , lik$outcome , "_missingness[u]] <- " , lik$outcome , "_impute[u]" , ";\n" )
                    else
                        m_tpars2 <- concat( m_tpars2 , indent , outcome , "[" , lik$outcome , "_missingness] <- " , lik$outcome , "_impute" , ";\n" )
                    # add x_impute to start list
                    imputer_name <- concat(lik$outcome,"_impute")
                    start[[ imputer_name ]] <- rep( impute_bank[[lik$outcome]]$init , impute_bank[[lik$outcome]]$N_miss )
                    if ( N_miss>1 )
                        types[[ imputer_name ]] <- c("vector", concat("[",lik$N_name,"_missing]") )
                    else
                        # only one missing value, so cannot use vector
                        types[[ imputer_name ]] <- "real"
                    # make sure x_impute is in pars
                    #pars[[ imputer_name ]] <- imputer_name
                }
                
                # ordinary sampling statement
                m_model_txt <- concat( m_model_txt , indent , outcome , " ~ " , lik$likelihood , "( " , parstxt , " )" , lik$T_text , ";\n" )
                
                # don't add deviance calc when imputed predictor
                if ( !(lik$outcome %in% names(impute_bank)) ) {
                    # get stan suffix
                    the_suffix <- templates[[lik$template]]$stan_suffix
                    m_gq <- concat( m_gq , indent , "dev <- dev + (-2)*" , lik$likelihood , the_suffix , "( " , outcome , " | " , parstxt , " )" , lik$T_text , ";\n" )
                }
                
            }
            
            # add N variable to data block
            the_N_name <- lik$N_name
            if ( i > 0 ) {
                if ( lik$outcome %in% names(impute_bank) ) {
                    the_N_name <- concat( the_N_name , "_missing" )
                    if ( i==1 ) {
                        warning(concat("Missing values in outcome variable ",lik$outcome," being imputed (predicted) through model. The model should fit fine. But helper functions like WAIC may not work."))
                    }
                }
                #m_data <- concat( m_data , indent , "int<lower=1> " , the_N_name , ";\n" )
                fp[['used_predictors']][[the_N_name]] <- list( var=the_N_name , type="index" )
                # warn about imputation on top-level outcome
                
            }
            
            # add number of cases to data list
            if ( lik$outcome %in% names(impute_bank) ) {
                d[[ the_N_name ]] <- as.integer( impute_bank[[ lik$outcome ]]$N_miss )
            } else {
                d[[ the_N_name ]] <- as.integer( lik$N_cases )
            }
            
        } # likelihood
        
    } # loop over fp_order
    
    # add number of cases to data list, in case no likelihood found
    d[[ "N" ]] <- as.integer( length(d[[1]]) )
    
    # compose generated quantities
    m_gq <- concat( m_model_declare , indent , "real dev;\n" , indent , "dev <- 0;\n" , m_gq )
    
    # general data length data declare
    #m_data <- concat( m_data , indent , "int<lower=1> " , "N" , ";\n" )
    fp[['used_predictors']][['N']] <- list( var="N" , type="index" )
    # data from used_predictors list
    n <- length( fp[['used_predictors']] )
    if ( n > 0 ) {
        for ( i in 1:n ) {
            var <- fp[['used_predictors']][[i]]
            type <- "real"
            # integer check
            if ( class(d[[var$var]])=="integer" ) type <- "int"
            if ( class(d[[var$var]])=="matrix" ) {
                type <- "matrix"
                if ( is.null(var$N) ) var$N <- dim(d[[var$var]])
            }
            # coerce outcome type
            if ( !is.null(var$type) ) type <- var$type
            # build
            if ( type=="matrix" ) {
                # rely on var$N being vector of matrix dimensions
                m_data <- concat( m_data , indent , type , "[",var$N[1],",",var$N[2],"] " , var$var , ";\n" )
            }
            if ( type=="index" ) {
                # index variable go at top of m_data
                m_data <- concat( indent , "int<lower=1> " , var$var , ";\n" , m_data )
            }  
            if ( type!="index" & type!="matrix" ) {
                m_data <- concat( m_data , indent , type , " " , var$var , "[" , var$N , "];\n" )
            }
        }#i
    }
    
    # declare parameters from start list
    # use any custom constraints in constraints list
    # first merge passed start with start built from priors
    if ( length(start_prior[[1]])>0 ) {
        start_p2 <- vector(mode="list",length=chains)
        for( ch in 1:chains ) start_p2[[ch]] <- list()
        n <- length(start_prior[[1]])
        # reverse index order, so parameters appear in same order as formula
        # need to do this, as we passed back-to-front when parsing formula
        for ( ki in n:1 ) {
            parname <- names(start_prior[[1]])[ki]
            for ( ch in 1:chains ) {
                start_p2[[ch]][[parname]] <- start_prior[[ch]][[parname]]
            }#ch
        } #ki
        # now merge start with start_p2
        # need to clone start first, then use unlist to merge each chain
        # start_p2 (inits sampled from priors) is unique to each chain at this point
        start_cloned <- vector(mode="list",length=chains)
        for ( ch in 1:chains ) start_cloned[[ch]] <- start
        start <- start_cloned
        for ( ch in 1:chains )
            start[[ch]] <- unlist( list( start_cloned[[ch]] , start_p2[[ch]] ) , recursive=FALSE )
    } else {
        # no inits sampled from priors
        # so just build the explicit start list
        start_cloned <- vector(mode="list",length=chains)
        for ( ch in 1:chains ) start_cloned[[ch]] <- start
        start <- start_cloned
    }

    ################ ATTENTION ###################
    # from this point on, start is a list of lists

    n <- length( start[[1]] )
    if ( n > 0 ) {
        for ( i in 1:n ) {
            pname <- names(start[[1]])[i]
            type <- "real"
            type_dim <- ""
            constraint <- ""
            
            if ( class(start[[1]][[i]])=="matrix" ) {
                # check for square matrix? just use nrow for now.
                #type <- concat( "cov_matrix[" , nrow(start[[i]]) , "]" )
                type <- "cov_matrix"
                type_dim <- concat( "[" , nrow(start[[1]][[i]]) , "]" )
                # corr_matrix check by naming convention
                #Rho_check <- grep( "Rho" , pname )
                #if ( length(Rho_check)>0 ) type <- concat( "corr_matrix[" , nrow(start[[i]]) , "]" )
            }
            
            # add correct length to numeric vectors (non-matrix)
            if ( type=="real" & length(start[[1]][[i]])>1 ) {
                #type <- concat( "vector[" , length(start[[i]]) , "]" )
                type <- "vector"
                type_dim <- concat( "[" , length(start[[1]][[i]]) , "]" )
                #if ( length(grep("sigma",pname))>0 )
                    #type <- concat( "vector<lower=0>[" , length(start[[i]]) , "]" )
                # check for varying effect vector by peeking at vprior list
                # we want symbolic length name here, if varying effect vector
                nvp <- length( fp[['vprior']] )
                if ( nvp > 0 ) {
                    for ( j in 1:nvp ) {
                        pars_out <- fp[['vprior']][[j]]$pars_out
                        for ( k in 1:length(pars_out) ) {
                            if ( pars_out[[k]]==pname ) {
                                #type <- concat( "vector[N_" , fp[['vprior']][[j]]$group , "]" )
                                type <- "vector"
                                type_dim <- concat( "[N_" , fp[['vprior']][[j]]$group , "]" )
                            }
                        }#k
                    }#j
                }
            }
            
            # add non-negative restriction to any parameter with 'sigma' in name
            #if ( length(grep("sigma",pname))>0 ) {
            #    if ( type=="real" ) constraint <- "<lower=0>"
            #}
            
            # any custom constraint?
            constrainttxt <- constraints[[pname]]
            if ( !is.null(constrainttxt) ) {
                constraint <- concat( "<" , constrainttxt , ">" )
                # check for positive constrain and validate start value
                # this is needed in case implicitly using half-distributions
                if ( constrainttxt == "lower=0" ) {
                    for ( ch in 1:chains )
                        start[[ch]][[pname]] <- abs( start[[ch]][[pname]] )
                }
            }
            
            # any custom type?
            mytype <- types[[pname]]
            if ( !is.null(mytype) ) {
                if ( length(mytype)==1 )
                    type <- mytype
                else {
                    # must have custom dims in [2]
                    type <- mytype[1]
                    type_dim <- mytype[2]
                }
            }
            
            # add to parameters block
            m_pars <- concat( m_pars , indent , type , constraint , type_dim , " " , pname , ";\n" )
        }#i
    }
    
    # compose any transformed parameters in transpars
    if ( length(transpars) > 0 ) {
        for ( i in 1:length(transpars) ) {
            if ( length(transpars[[i]])>1 ) {
                # declaration [1] and code [2]
                tptxt <- concat( indent , transpars[[i]][1] , ";\n" )
                m_tpars1 <- concat( m_tpars1 , tptxt )
                tptxt <- concat( indent , transpars[[i]][2] , ";\n" )
                m_tpars2 <- concat( m_tpars2 , tptxt )
            } else {
                # just code? --- not a used case yet
            }
        }#i
    }# transpars > 0
    
    # put it all together
    
    # function to add header/footer to each block in code
    # empty blocks remain empty
    blockify <- function(x,header,footer) {
        if ( x != "" ) x <- concat( header , x , footer )
        return(x)
    }
    m_data <- blockify( m_data , "data{\n" , "}\n" )
    m_pars <- blockify( m_pars , "parameters{\n" , "}\n" )
    m_tpars1 <- blockify( m_tpars1 , "transformed parameters{\n" , "" )
    m_tpars2 <- blockify( m_tpars2 , "" , "}\n" )
    m_gq <- blockify( m_gq , "generated quantities{\n" , "}\n" )
    
    #model_code <- concat( m_data , m_pars , m_tpars1 , m_tpars2 , "model{\n" ,  m_model_declare , m_model_priors , m_model_lm , m_model_lik , "}\n" , m_gq )
    
    model_code <- concat( m_data , m_pars , m_tpars1 , m_tpars2 , "model{\n" ,  m_model_declare , m_model_txt , "}\n" , m_gq )

    # from rstan version 2.10.0: change all '<-' to '='
    model_code <- gsub( "<-" , "=" , model_code , fixed=TRUE )
    
    if ( debug==TRUE ) cat(model_code)

##############################
# end of Stan code compilation
##############################
    
    ########################################
    # fit model
    
    # build pars vector
    # use ours, unless user provided one
    if ( missing(pars) ) {
        pars <- names(start[[1]])
        elected <- names(pars_elect)
        pars <- c( pars , elected , "dev" )
        if ( length(pars_hide)>0 ) {
            exclude_idx <- which( pars %in% names(pars_hide) )
            pars <- pars[ -exclude_idx ]
        }
    }
    
    if ( sample==TRUE ) {
        #require(rstan) # don't need anymore
        
        # sample
        modname <- deparse( flist[[1]] )
        #initlist <- list()
        #for ( achain in 1:chains ) initlist[[achain]] <- start
        initlist <- start # should have one list for each chain

        if ( missing(rng_seed) ) rng_seed <- sample( 1:1e5 , 1 )
        if ( is.null(previous_stanfit) )
            fit <- try(stan( model_code=model_code , model_name=modname , data=d , init=initlist , iter=iter , warmup=warmup , chains=chains , cores=cores , pars=pars , seed=rng_seed , ... ))
        else
            fit <- try(stan( fit=previous_stanfit , model_name=modname , data=d , init=initlist , iter=iter , warmup=warmup , chains=chains , cores=cores , pars=pars , seed=rng_seed , ... ))

        if ( class(fit)=="try-error" ) {
                # something went wrong in at least one chain
                msg <- attr(fit,"condition")$message
                if ( cores > 1 )
                    stop(concat("Something went wrong in at least one chain. Debug your model while setting chains=1 and cores=1. Once the model is working with a single chain and core, try using multiple chains/cores again.\n",msg))
                else
                    stop(concat("Something went wrong, when calling Stan. Check any debug messages for clues, detective.\n",msg))
            }
        
        ###### PREVIOUS multicore code, before Stan added its own ########
        if ( FALSE ) {
        if ( chains==1 | cores==1 ) {
        
            # single core
            # so just run the model
            if ( missing(rng_seed) ) rng_seed <- sample( 1:1e5 , 1 )
            if ( is.null(previous_stanfit) )
                fit <- stan( model_code=model_code , model_name=modname , data=d , init=initlist , iter=iter , warmup=warmup , chains=chains , pars=pars , seed=rng_seed , ... )
            else
                fit <- stan( fit=previous_stanfit , model_name=modname , data=d , init=initlist , iter=iter , warmup=warmup , chains=chains , pars=pars , seed=rng_seed , ... )
            
        } else {
        
            # multi-core, multi-chain
            # so pre-compile then run again
            # note use of chain number 1 inits only here
            # not sure inits are used at all, given chains=0
            if ( is.null(previous_stanfit) ) {
                message("Precompiling model...")
                fit_prep <- stan( model_code=model_code , model_name=modname , data=d , init=start[[1]] , iter=1 , chains=0 , pars=pars , test_grad=FALSE , refresh=-1 , ... )
            } else {
                # reuse compiled model
                fit_prep <- previous_stanfit
            }
            
            # if no seed provided, make one
            if ( missing(rng_seed) ) rng_seed <- sample( 1:1e5 , 1 )
            
            message(concat("Dispatching ",chains," chains to ",cores," cores."))
            sys <- .Platform$OS.type
            if ( sys=='unix' ) {
                # Mac or Linux
                # hand off to mclapply
                sflist <- mclapply( 1:chains , mc.cores=cores ,
                    function(chainid)
                        stan( fit=fit_prep , data=d , init=list(start[[chainid]]) , pars=pars , iter=iter , warmup=warmup , chains=1 , seed=rng_seed , chain_id=chainid , ... )
                )
            } else {
                # Windows
                # so use parLapply instead
                CL = makeCluster(cores)
                #data <- object@data
                env0 <- list( fit_prep=fit_prep, d=d, pars=pars, rng_seed=rng_seed, iter=iter, warmup=warmup , start=start )
                clusterExport(cl = CL, 
                    c("fit_prep","d","pars","iter","warmup","rng_seed","start") , as.environment(env0) )
                sflist <- parLapply(CL, 1:chains, fun = function(cid) {
					require(rstan) # will CRAN tolerate this?
                    stan( fit=fit_prep , data = d, pars = pars, chains = 1, 
                      iter = iter, warmup = warmup, seed = rng_seed, 
                      chain_id = cid, init=list(start[[cid]]) )
                })
            }
            # merge result
            fit <- try( sflist2stanfit(sflist) )
            if ( class(fit)=="try-error" ) {
                # something went wrong in at least one chain
                stop("Something went wrong in at least one chain. Debug your model while setting chains=1. Once the model is working with a single chain, try using multiple chains again.")
            }
            
        }#multicore
        }
        # END PREVIOUS multicore code
        #############################
        
    } else {
        fit <- NULL
    }
    
    ########################################
    # build result
    
    coef <- NULL
    varcov <- NULL
    fp$impute_bank <- impute_bank
    if ( sample==TRUE ) {
        # compute expected values of parameters
        s <- summary(fit)$summary
        s <- s[ -which( rownames(s)=="lp__" ) , ]
        s <- s[ -which( rownames(s)=="dev" ) , ]
        if ( !is.null(dim(s)) ) {
            coef <- s[,1]
            # compute variance-covariance matrix
            varcov <- matrix(NA,nrow=nrow(s),ncol=nrow(s))
            diag(varcov) <- s[,3]^2
        } else {
            coef <- s[1]
            varcov <- matrix( s[3]^2 , 1 , 1 )
            names(coef) <- names(start[[1]])
        }
        
        # compute DIC
        dev.post <- extract(fit, "dev", permuted = TRUE, inc_warmup = FALSE)
        dbar <- mean( dev.post$dev )
        # to compute dhat, need to feed parameter averages back into compiled stan model
        post <- extract( fit )
        Epost <- list()
        for ( i in 1:length(post) ) {
            dims <- length( dim( post[[i]] ) )
            name <- names(post)[i]
            if ( name!="lp__" & name!="dev" ) {
                if ( dims==1 ) {
                    Epost[[ name ]] <- mean( post[[i]] )
                } else {
                    Epost[[ name ]] <- apply( post[[i]] , 2:dims , mean )
                }
            }
        }#i
        
        if ( debug==TRUE ) print( Epost )
        
        # push expected values back through model and fetch deviance
        #message("Taking one more sample now, at expected values of parameters, in order to compute DIC")
        #fit2 <- stan( fit=fit , init=list(Epost) , data=d , pars="dev" , chains=1 , iter=2 , refresh=-1 , cores=1 )
        fit2 <- sampling( fit@stanmodel , init=list(Epost) , data=d , pars="dev" , chains=1 , iter=1 , cores=1 )
        dhat <- as.numeric( extract(fit2,"dev") )
        pD <- dbar - dhat
        dic <- dbar + pD
        
        # if (debug==TRUE) print(Epost)
        
        # build result
        result <- new( "map2stan" , 
            call = match.call(), 
            model = model_code,
            stanfit = fit,
            coef = coef,
            vcov = varcov,
            data = d,
            start = start,
            pars = pars,
            formula = flist.orig,
            formula_parsed = fp )
        
        attr(result,"df") = length(result@coef)
        attr(result,"DIC") = dic
        attr(result,"pD") = pD
        attr(result,"deviance") = dhat
        try( 
            if (!missing(d)) attr(result,"nobs") = length(d[[ fp[['likelihood']][[1]][['outcome']] ]]) , 
            silent=TRUE
        )
        
        # compute WAIC?
        if ( WAIC==TRUE ) {
            message("Computing WAIC")
            waic <- try(WAIC( result , n=0 , pointwise=TRUE )) # n=0 to use all available samples
            attr(result,"WAIC") = waic
        }
        
        # check divergent iterations
        nd <- divergent(fit)
        if ( nd > 0 ) {
            warning( concat("There were ",nd," divergent iterations during sampling.\nCheck the chains (trace plots, n_eff, Rhat) carefully to ensure they are valid.") )
        }
        
    } else {
        # just return list
        result <- list(
            call = match.call(), 
            model = model_code,
            data = d,
            start = start,
            pars = pars,
            formula = flist.orig,
            formula_parsed = fp )
    }
    
    return( result )
    
}

# EXAMPLES
if ( FALSE ) {

library(rethinking)

# simulate data
library(MASS)
N <- 500 # 1000 cases
J <- 20 # 100 clusters
J2 <- 10
NperJ <- N/J
sigma <- 2 # top-level standard deviation
mu <- c(10,-0.5) # means of varying effects coefficients
x <- runif(N,min=-2,max=2) # predictor variable
x2 <- runif(N,min=-2,max=2)
id <- rep( 1:J , each=NperJ ) # cluster id's
id2 <- rep( 1:J2 , each=N/J2 )
Sigma <- matrix( 0 , nrow=2 , ncol=2 ) # var-cov matrix
Sigma[1,1] <- 2
Sigma[2,2] <- 0.2
Sigma[1,2] <- Sigma[2,1] <- -0.8 * sqrt( Sigma[1,1] * Sigma[2,2] )
beta <- mvrnorm( J , mu=mu , Sigma=Sigma )
y <- rnorm( N , mean=beta[id,1]+beta[id,2]*x , sd=sigma )
y2 <- rbinom( N , size=1 , prob=logistic( y-8 ) )

# fitting tests


# cross classified
f <- list(
    y ~ dnorm(mu,sigma),
    mu ~ a + aj1 + aj2 + b*x,
    aj1|id ~ dnorm( 0 , sigma_id ),
    aj2|id2 ~ dnorm( 0 , sigma_id2 ),
    a ~ dnorm(0,10),
    b ~ dnorm(0,1),
    sigma ~ dcauchy(0,1),
    sigma_id ~ dcauchy(0,1),
    sigma_id2 ~ dcauchy(0,1)
)
startlist <- list(a=10,b=0,sigma=2,sigma_id=1,sigma_id2=1,aj1=rep(0,J),aj2=rep(0,J2))
m <- map2stan( f , data=list(y=y,x=x,id=id,id2=id2) , start=startlist , sample=TRUE , debug=FALSE )


# random slopes with means inside multi_normal
f <- list(
    y ~ dnorm(mu,sigma),
    mu ~ aj + bj*x,
    c(aj,bj)|id ~ dmvnorm( c(a,b) , Sigma_id ),
    a ~ dnorm(0,10),
    b ~ dnorm(0,1),
    sigma ~ dcauchy(0,1),
    Sigma_id ~ inv_wishart(3,diag(2)) # or dinvwishart
)
startlist <- list(a=10,b=0,sigma=2,Sigma_id=diag(2),aj=rep(0,J),bj=rep(0,J))
m2 <- map2stan( f , data=list(y=y,x=x,id=id) , start=startlist , sample=TRUE , debug=FALSE )

cat(m$model)


# 
f2 <- list(
    y ~ dnorm(mu,sigma),
    mu ~ a + aj + b*x,
    aj|id ~ dnorm( 0 , sigma_a ),
    a ~ dnorm(0,10),
    b ~ dnorm(0,1),
    sigma ~ dcauchy(0,1),
    sigma_a ~ dcauchy(0,1)
)
startlist <- list(a=10,b=0,sigma=3,sigma_a=1,aj=rep(0,J))
m <- map2stan( f2 , data=list(y=y,x=x,id=id) , start=startlist , sample=TRUE , debug=FALSE )

# now with fixed effect inside prior
f4 <- list(
    y ~ dnorm(mu,sigma),
    mu ~ aj + b*x,
    aj|id ~ dnorm( a , sigma_a ),
    a ~ dnorm(0,10),
    b ~ dnorm(0,1),
    sigma ~ dcauchy(0,1),
    sigma_a ~ dcauchy(0,1)
)
startlist <- list(a=10,b=0,sigma=2,sigma_a=1,aj=rep(0,J))
m2 <- map2stan( f4 , data=list(y=y,x=x,id=id) , start=startlist , sample=TRUE , debug=FALSE )

# random slopes
f <- list(
    y ~ dnorm(mu,sigma),
    mu ~ a + aj + (b+bj)*x,
    c(aj,bj)|id ~ dmvnorm( 0 , Sigma_id ),
    a ~ dnorm(0,10),
    b ~ dnorm(0,1),
    sigma ~ dcauchy(0,1),
    Sigma_id ~ inv_wishart(3,diag(2)) # or dinvwishart
)
startlist <- list(a=10,b=0,sigma=2,Sigma_id=diag(2),aj=rep(0,J),bj=rep(0,J))
m <- map2stan( f , data=list(y=y,x=x,id=id) , start=startlist , sample=TRUE , debug=FALSE )


f3 <- list(
    y ~ dbinom(1,theta),
    logit(theta) ~ a + aj + b*x,
    aj|id ~ dnorm( 0 , sigma_a ),
    a ~ dnorm(0,10),
    b ~ dnorm(0,1),
    sigma_a ~ dcauchy(0,1)
)



f5 <- list(
    y ~ dnorm(mu,sigma),
    mu ~ aj + bj*x,
    c(aj,bj)|id ~ dmvnorm( c(a,b) , Sigma_id ),
    a ~ dnorm(0,10),
    b ~ dnorm(0,1),
    sigma ~ dcauchy(0,1)
)


} #EXAMPLES
