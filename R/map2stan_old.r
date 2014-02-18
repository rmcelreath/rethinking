# map2stan - convert map-like formula to Stan model code
# TO-DO:
# (-) gamma_poisson parameter conversion?
# (-) allow more than one likelihood/linear model --- need to mark N variables
# (-) check for NAs in data
# (-) Allow explicit start values for sigma's in varying effects
# (-) Allow optional random start values for all parameters
# (-) fix plotchains() error with random slopes
# (-) allow covmatrix priors on standard devs and corrmatrix (sigma*Rho*sigma?)
#     Could just allow specifying dmvnorm( 0 , diag(sigma)%*%Rho%*%diag(sigma) )
#     Then user can put priors on sigma vector and Rho matrix

map2stan_old <- function( flist , data , start , sample=TRUE , iter=2000 , chains=1 , debug=FALSE , ... ) {
    
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
    if ( missing(data) ) stop( "'data' required." )
    if ( !( class(data) %in% c("list","data.frame") ) ) {
        stop( "'data' must be of class list or data.frame." )
    }
    
    flist.orig <- flist
    start.orig <- start
    
    ########################################
    # check for common issues
    # (1) NAs in data vectors
    # tricky, because only want to check used variables
    
    # (2) type mismatch for outcome variable
    
    ########################################
    # private functions
    
    # for converting characters not allows by Stan
    undot <- function( astring ) {
        astring <- gsub( "." , "_" , astring , fixed=TRUE )
        astring
    }
    
    # substitution list
    density_sub_list <- list(
        dnorm = 'normal',
        normal = 'normal',
        dmvnorm = 'multi_normal',
        multi_normal = 'multi_normal',
        dunif = 'uniform',
        uniform = 'uniform',
        dcauchy = 'cauchy',
        cauchy = 'cauchy',
        dexp = 'exp',
        exp = 'exp',
        dbinom = 'binomial',
        binomial = 'binomial',
        dpois = 'poisson',
        poisson = 'poisson',
        dgamma = 'gamma',
        gamma = 'gamma',
        dlaplace = 'double_exponential',
        double_exponential = 'double_exponential',
        laplace = 'double_exponential',
        dbetabinom = 'beta_binomial',
        beta_binomial = 'beta_binomial'
    )
    
    inverse_links <- list(
        log = 'exp',
        logit = 'inv_logit',
        logodds = 'inv_logit'
    )
    
    ########################################
    # parse formulas
    
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
    
    # extract likelihood(s)
    extract_likelihood <- function( fl ) {
        outcome <- as.character( fl[[2]] )
        likelihood <- fl[[3]][[1]]
        likelihood <- density_sub_list[[ as.character(likelihood) ]]
        likelihood_pars <- list()
        for ( i in 1:(length(fl[[3]])-1) ) likelihood_pars[[i]] <- as.character( fl[[3]][[i+1]] )
        # result
        list( 
            outcome = undot(outcome) ,
            likelihood = likelihood ,
            pars = likelihood_pars
        )
    }
    
    # extract linear model(s)
    extract_linearmodel <- function( fl ) {
        # check for link function
        if ( length(fl[[2]])==1 ) {
            parameter <- as.character( fl[[2]] )
            link <- "identity"
        } else {
            # link!
            parameter <- as.character( fl[[2]][[2]] )
            link <- as.character( fl[[2]][[1]] )
        }
        RHS <- paste( deparse(fl[[3]]) , collapse=" " )
        list(
            parameter = parameter ,
            RHS = RHS ,
            index = 'i' ,
            link = link
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
        
        density <- density_sub_list[[ as.character(fl[[3]][[1]]) ]]
        # must be either normal or multi_normal, so always 2 parameters
        mu <- fl[[3]][[2]]
        Sigma <- as.character( fl[[3]][[3]] )
        pars_in <- list( mu , Sigma )
        list(
            pars_out = pars ,
            density = density ,
            pars_in = pars_in ,
            group = group_var
        )
    }
    
    # extract ordinary prior
    extract_prior <- function( fl ) {
        apar <- as.character( fl[[2]] )
        adensity <- density_sub_list[[ as.character(fl[[3]][[1]]) ]]
        inpars <- list()
        n <- length( fl[[3]] ) - 1
        for ( i in 1:n ) {
            inpars[[i]] <- fl[[3]][[i+1]]
        }
        list(
            par_out = apar ,
            density = adensity ,
            pars_in = inpars
        )
    }
    
    # build parsed list
    fp <- list( 
        likelihood = list() ,
        lm = list() ,
        vprior = list(),
        prior = list()
    )
    fp[['used_predictors']] <- list()
    d <- as.list(data)
    
    # pass over formulas and extract into correct slots
    for ( i in 1:length(flist) ) {
        # test for likelihood
        if ( i==1 ) {
            lik <- extract_likelihood( flist[[i]] )
            fp[['likelihood']][[1]] <- lik
            # check for binomial size variable and mark used
            if ( lik$likelihood=='binomial' | lik$likelihood=='beta_binomial' ) {
                sizename <- lik$pars[[1]]
                if ( !is.null( d[[sizename]] ) ) {
                    nup <- length(fp[['used_predictors']])
                    fp[['used_predictors']][[nup+1]] <- undot(sizename)
                }
            }
            next
        }
        # test for linear model
        is_linearmodel <- FALSE
        if ( length(flist[[i]][[3]]) > 1 ) {
            fname <- as.character( flist[[i]][[3]][[1]] )
            if ( fname=="+" | fname=="*" | fname=="-" | fname=="/" | fname=="%*%" | fname %in% c('exp','logistic') ) {
                is_linearmodel <- TRUE
            }
        } else {
            # RHS is length 1, so can't be a prior or density statement
            # must be a linear model with a single symbol in it, like "p ~ a"
            is_linearmodel <- TRUE
        }
        if ( is_linearmodel==TRUE ) {
            n <- length( fp[['lm']] )
            fp[['lm']][[n+1]] <- extract_linearmodel( flist[[i]] )
            next
        }
        
        # test for varying effects prior
        if ( length( flist[[i]][[2]] ) == 3 ) {
            fname <- as.character( flist[[i]][[2]][[1]] )
            if ( fname=="|" ) {
                n <- length( fp[['vprior']] )
                fp[['vprior']][[n+1]] <- extract_vprior( flist[[i]] )
                next
            }
        }
        # an ordinary prior?
        n <- length( fp[['prior']] )
        fp[['prior']][[n+1]] <- extract_prior( flist[[i]] )
    }
    
    # add index brackets in linear models
    
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
    detectvar <- function( target , x ) {
        wild <- "[()=*+ ]"
        x2 <- paste( " " , x , " " , collapse="" , sep="" )
        pattern <- paste( wild , target , wild , sep="" , collapse="" )
        m <- regexpr( pattern , x2 )
        return(m)
    }
    
    # go through linear models
    for ( i in 1:length(fp[['lm']]) ) {
        # for each variable in data, add index
        index <- fp[['lm']][[i]][['index']]
        for ( v in names(d) ) {
            # tag if used
            used <- detectvar( v , fp[['lm']][[i]][['RHS']] )
            if ( used > -1 ) {
                nup <- length(fp[['used_predictors']])
                fp[['used_predictors']][[nup+1]] <- undot(v)
            }
            # add index and undot the name
            fp[['lm']][[i]][['RHS']] <- indicize( v , index , fp[['lm']][[i]][['RHS']] , replace=undot(v) )
        }
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
                        jindexn <- paste( jindex , "," , k , collapse="" , sep="" ) 
                        fp[['lm']][[i]][['RHS']] <- indicize( var , jindexn , fp[['lm']][[i]][['RHS']] , vname )
                    }
                }
            }
        }
    }
    
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
    concat <- function( ... ) {
        paste( ... , collapse="" , sep="" )
    }
    
    # data block
    
    m_data <- "data{\n"
    m_data <- concat( m_data , indent , "int N;\n" )
    # add outcome(s)
    for ( i in 1:length(fp[['likelihood']]) ) {
        oname <- fp[['likelihood']][[i]][['outcome']]
        dtype <- ifelse( class(d[[oname]])=="numeric" , "real" , "int" )
        m_data <- concat( m_data , indent , dtype , " " , oname , "[N];\n" )
        d[["N"]] <- length( d[[oname]] )
    }
    # add predictors
    if ( debug==TRUE ) print(fp[['used_predictors']])
    np <- length(fp[['used_predictors']])
    if ( np > 0 ) {
        for ( i in 1:np ) {
            vname <- fp[['used_predictors']][[i]]
            dtype <- ifelse( class(d[[vname]])=="numeric" , "real" , "int" )
            m_data <- concat( m_data , indent , dtype , " " , vname , "[N];\n" )
        } #i
    } #np
    # add index variables and group sample sizes
    nv <- length(fp[['vprior']])
    if ( nv > 0 ) {
        for ( i in 1:nv ) {
            gname <- fp[['vprior']][[i]][['group']]
            Nname <- concat( "N_" , gname )
            m_data <- concat( m_data , indent , "int" , " " , gname , "[N];\n" )
            m_data <- concat( m_data , indent , "int" , " " , Nname , ";\n" )
            # add sample size to data list
            d[[Nname]] <- length( unique( d[[gname]] ) )
        }
    } #nv
    m_data <- concat( m_data , "}\n" )
    
    # parameters block
    
    # use start list to build in fixed effects
    m_pars <- "parameters{\n"
    for ( i in 1:length(start) ) {
        pname <- names(start)[i]
        dtype <- "real"
        if ( pname=="sigma" ) dtype <- "real<lower=0>"
        if ( pname=="theta" & fp[['likelihood']][[1]][['likelihood']]=="beta_binomial" ) dtype <- "real<lower=0>"
        m_pars <- concat( m_pars , indent , dtype , " " , pname , ";\n" )
    }
    # add varying effects
    if ( nv > 0 ) {
        for ( i in 1:nv ) {
            gname <- fp[['vprior']][[i]][['group']]
            Nname <- concat( "N_" , gname )
            num_pars <- length( fp[['vprior']][[i]][['pars_out']] )
            if ( num_pars > 1 ) {
                m_pars <- concat( m_pars , indent , "vector[" , num_pars , "] " , "vary_" , gname , "[" , Nname , "];\n" )
                Sigmaname <- fp[['vprior']][[i]][['pars_in']][[2]]
                m_pars <- concat( m_pars , indent , "cov_matrix[" , num_pars , "] " , Sigmaname , ";\n" )
            } else {
                outname <- fp[['vprior']][[i]][['pars_out']][[1]]
                m_pars <- concat( m_pars , indent , "real " , outname , "[" , Nname , "];\n" )
                sigmaname <- fp[['vprior']][[i]][['pars_in']][[2]]
                m_pars <- concat( m_pars , indent , "real<lower=0> " , sigmaname , ";\n" )
            }
        }#i
    }#nv
    m_pars <- concat( m_pars , "}\n" )
    
    # transformed parameters
    # these are used just to translate matrix vary_ parameters back to the vector parameters in the user formula
    nv <- length( fp[['vprior']] )
    m_tpar <- "transformed parameters{\n"
    use_tpar <- FALSE
    if ( nv > 0 ) {
        for ( i in 1:nv ) {
            pars <- fp[['vprior']][[i]][['pars_out']]
            parsin <- fp[['vprior']][[i]][['pars_in']]
            if ( length(pars)>1 ) {
                use_tpar <- TRUE
                # declare vector outputs
                for ( j in 1:length(pars) ) {
                    m_tpar <- concat( m_tpar , indent , "vector[" , Nname , "] " , pars[[j]] , ";\n" )
                } #j
                # check for parameters in first input -- may need transform there
                use_input_transform <- FALSE
                if ( fp[['vprior']][[i]][['density']]=="multi_normal" ) {
                    if ( length(parsin[[1]]) > 1 ) {
                        use_input_transform <- TRUE
                        # vector of parameters like c(a,b)
                        # will use transform mu_GROUP in model block
                        # so declare vector here
                        m_tpar <- concat( m_tpar , indent , "vector[" , length(pars) , "] mu_" , fp[['vprior']][[i]][['group']] , ";\n" )
                    }
                }
                # assign 'real' outputs to vector
                m_tpar <- concat( m_tpar , indent , "for ( j in 1:" , Nname , " ) {\n" )
                for ( j in 1:length(pars) ) {
                    fakename <- concat( "vary_" , fp[['vprior']][[i]][['group']] )
                    m_tpar <- concat( m_tpar , indent , indent , pars[[j]] , "[j] <- " , fakename , "[j," , j , "];\n" )
                } #j
                m_tpar <- concat( m_tpar , indent , "}\n" )
                if ( use_input_transform==TRUE ) {
                    # transform for vector of multi_normal means
                    mu <- parsin[[1]]
                    mu[[1]] <- NULL # erase function 'c'
                    vecname <- concat( "mu_" , fp[['vprior']][[i]][['group']] )
                    for ( j in 1:length(mu) ) {
                        m_tpar <- concat( m_tpar , indent , vecname , "[" , j , "] <- " , as.character(mu[[j]]) , ";\n" )
                    }
                }
            }
        }#nv
    }
    m_tpar <- concat( m_tpar , "}\n" )
    
    # model block and generated quantities block
    m_model <- "model{\n"
    m_gq <- concat( "generated quantities{\n" , indent , "real dev;\n" )
    
    # temporary variables --- names of linear models
    nf <- length( fp[['lm']] )
    if ( nf > 0 ) {
        for ( i in 1:nf ) {
            pname <- fp[['lm']][[i]][['parameter']]
            m_model <- concat( m_model , indent , "vector[N] " , pname , ";\n" )
            m_gq <- concat( m_gq , indent , "vector[N] " , pname , ";\n" )
        }
    } #nf
    m_gq <- concat( m_gq , indent , "dev <- 0;\n" )
    
    # priors
    np <- length( fp[['prior']] )
    if ( np > 0 ) {
        for ( i in 1:np ) {
            pname <- fp[['prior']][[i]][['par_out']]
            dname <- fp[['prior']][[i]][['density']]
            m_model <- concat( m_model , indent , pname , " ~ " , dname , "(" )
            np2 <- length(fp[['prior']][[i]][['pars_in']])
            m_model <- concat( m_model , " " , fp[['prior']][[i]][['pars_in']][[1]] )
            if ( np2 > 1 ) {
                for ( j in 2:np2 ) {
                    m_model <- concat( m_model , " , " , fp[['prior']][[i]][['pars_in']][[j]] )
                }#j
            }
            m_model <- concat( m_model , " );\n" )
        }#i
    } #np
    
    # varying effects
    nv <- length( fp[['vprior']] )
    if ( nv > 0 ) {
        for ( i in 1:nv ) {
            gname <- fp[['vprior']][[i]][['group']]
            Nname <- concat( "N_" , gname )
            dname <- fp[['vprior']][[i]][['density']]
            num_out <- length( fp[['vprior']][[i]][['pars_out']] )
            if ( num_out > 1 ) {
                # multivariate distribution
                if ( dname=="multi_normal" ) {
                    if ( length(fp[['vprior']][[i]][['pars_in']][[1]])>1 ) {
                        # a vector of inputs c() to a single density parameter
                        # try using transformed parameter vector for it
                        parsin <- fp[['vprior']][[i]][['pars_in']]
                        parsin[[1]] <- concat( "mu_" , gname )
                        inputstxt <- paste( parsin , collapse=" , " , sep="" )
                    }
                    if ( fp[['vprior']][[i]][['pars_in']][[1]]=="0" ) {
                        # repeat zero mean for vector mean of multi_normal
                        parsin <- fp[['vprior']][[i]][['pars_in']]
                        parsin[[1]] <- concat( "rep_vector(0," , num_out , ")" )
                        inputstxt <- paste( parsin , collapse=" , " , sep="" )
                    }
                } else {
                    # just copy inputs as given
                    inputstxt <- paste( fp[['vprior']][[i]][['pars_in']] , collapse=" , " , sep="" )
                }
                m_model <- concat( m_model , indent , "for ( j in 1:" , Nname , " ) " , "vary_" , gname , "[j] ~ " , dname , "( " , inputstxt , " );\n" )
            } else {
                # single variate distribution
                poname <- fp[['vprior']][[i]][['pars_out']][[1]]
                npi <- length( fp[['vprior']][[i]][['pars_in']] )
                if ( npi > 1 ) {
                    parsintxt <- paste( fp[['vprior']][[i]][['pars_in']] , collapse=" , " , sep="" )
                } else {
                    parsintxt <- fp[['vprior']][[i]][['pars_in']][[1]]
                }
                m_model <- concat( m_model , indent , "for ( j in 1:" , Nname , " ) " , poname , "[j] ~ " , dname , "( " , parsintxt , " );\n" )
            }
        }
    } #nv
    
    # linear models and likelihoods
    for ( i in 1:length(fp[['lm']]) ) {
        m_model <- concat( m_model , indent , "for ( i in 1:N ) {\n" )
        m_gq <- concat( m_gq , indent , "for ( i in 1:N ) {\n" )
        lhs <- fp[['lm']][[i]][['parameter']]
        rhs <- fp[['lm']][[i]][['RHS']]
        link <- fp[['lm']][[i]][['link']]
        m_model <- concat( m_model , indent , indent , lhs , "[i] <- " , rhs , ";\n" )
        m_gq <- concat( m_gq , indent , indent , lhs , "[i] <- " , rhs , ";\n" )
        if ( link != "identity" ) {
            # add inverse link
            invlink <- inverse_links[[link]]
            m_model <- concat( m_model , indent , indent , lhs , "[i] <- " , invlink , "(" , lhs , "[i]);\n" )
            m_gq <- concat( m_gq , indent , indent , lhs , "[i] <- " , invlink , "(" , lhs , "[i]);\n" )
        }
        # close N loop
        m_model <- concat( m_model , indent , "}\n" )
        m_gq <- concat( m_gq , indent , "}\n" )
        # likelihood
        outcome <- fp[['likelihood']][[i]][['outcome']]
        dname <- fp[['likelihood']][[i]][['likelihood']]
        pars <- fp[['likelihood']][[i]][['pars']]
        # special pars transforms to map from R to Stan functions
        if ( dname=="beta_binomial" ) {
            # p,theta to alpha,beta
            name.p <- pars[[2]]
            name.theta <- pars[[3]]
            pars[[2]] <- concat( name.p , "*" , name.theta )
            pars[[3]] <- concat( "(1-" , name.p , ")*" , name.theta )
        }
        # format pars
        pars <- paste( pars , collapse=" , " )
        # use vectorized form
        m_model <- concat( m_model , indent , outcome , " ~ " , dname , "( " , pars , " );\n" )
        m_gq <- concat( m_gq , indent , "dev <- (-2)*" , dname , "_log( " , outcome , " , " , pars , " );\n" )
    }
    m_model <- concat( m_model , "}\n" )
    m_gq <- concat( m_gq , "}\n" )
    
    # put it all together
    if ( use_tpar==FALSE ) m_tpar <- ""
    model_code <- concat( m_data , m_pars , m_tpar , m_model , m_gq )
    if ( debug==TRUE ) cat(model_code)
    
    ########################################
    # fit model
    
    if ( sample==TRUE ) {
        require(rstan)
        
        # add varying effects to start list
        # also build pars vector
        pars <- names(start.orig)
        pars <- c( pars , "dev" )
        nv <- length(fp[['vprior']])
        if ( nv > 0 ) {
            for ( i in 1:nv ) {
                # out pars
                npars <- length( fp[['vprior']][[i]][['pars_out']] )
                gname <- fp[['vprior']][[i]][['group']]
                ngroups <- length( unique( d[[ gname ]] ) )
                if ( npars > 1 ) {
                    # vector of out parameters
                    parname <- concat( "vary_" , fp[['vprior']][[i]][['group']] )
                    start[[parname]] <- matrix(0,nrow=ngroups,ncol=npars)
                    for ( j in 1:npars ) {
                        apar <- fp[['vprior']][[i]][['pars_out']][[j]]
                        pars <- c( pars , apar )
                    }#j
                } else {
                    # single out parameter
                    parname <- fp[['vprior']][[i]][['pars_out']][[1]]
                    start[[parname]] <- rep(0,ngroups)
                    pars <- c( pars , parname )
                }
                # in pars
                Sname <- fp[['vprior']][[i]][['pars_in']][[2]]
                if ( npars > 1 )
                    start[[Sname]] <- diag(npars)
                else
                    start[[Sname]] <- 1
                pars <- c( pars , Sname )
            }#i
        }#nv
        
        if ( debug==TRUE ) print( start )
        
        # sample
        modname <- deparse( flist[[1]] )
        initlist <- list()
        for ( achain in 1:chains ) initlist[[achain]] <- start
        if ( use_tpar==FALSE ) {
            fit <- stan( model_code=model_code , model_name=modname , data=d , init=initlist , iter=iter , chains=chains , pars=pars , ... )
        } else {
            # use pars vector to omit matrix varying effect parameters
            fit <- stan( model_code=model_code , model_name=modname , data=d , init=initlist , iter=iter , chains=chains , pars=pars , ... )
        }
    } else {
        fit <- NULL
    }
    
    ########################################
    # build result
    
    coef <- NULL
    varcov <- NULL
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
            names(coef) <- names(start.orig)
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
        
        # need to swap out parameter names for internal vary_ names, if any
        nv <- length( fp[['vprior']] )
        if ( nv > 0 ) {
            for ( i in 1:nv ) {
                pout <- fp[['vprior']][[i]][['pars_out']]
                if ( length(pout) > 1 ) {
                    num_groups <- length( Epost[[ pout[[1]] ]] )
                    pnew <- matrix( 0 , nrow=num_groups , ncol=length(pout) )
                    for ( j in 1:length(pout) ) {
                        pnew[,j] <- Epost[[ pout[[j]] ]]
                        Epost[[ pout[[j]] ]] <- NULL
                    }#j
                    parname <- concat( "vary_" , fp[['vprior']][[i]][['group']] )
                    Epost[[ parname ]] <- pnew
                }
            }#i
        }#nv
        
        if ( debug==TRUE ) print( Epost )
        
        # push expected values back through model and fetch deviance
        message("Taking one more sample now, at expected values of parameters, in order to compute DIC:")
        fit2 <- stan( fit=fit , init=list(Epost) , data=d , pars="dev" , chains=1 , iter=1 )
        dhat <- as.numeric( extract(fit2,"dev") )
        pD <- dbar - dhat
        dic <- dbar + pD
        
        # if (debug==TRUE) print(Epost)
        
        # build result
        result <- new( "map2stan" , 
            call = match.call(), 
            stanfit = fit,
            coef = coef,
            vcov = varcov,
            data = d,
            start = start,
            formula = flist,
            formula_parsed = fp )
        attr(result,"df") = length(result@coef)
        attr(result,"DIC") = dic
        attr(result,"pD") = pD
        attr(result,"deviance") = dhat
        if (!missing(d)) attr(result,"nobs") = length(d[[ fp[['likelihood']][[1]][['outcome']] ]])
    } else {
        # just return list
        result <- list( 
            call = match.call(), 
            model = model_code,
            data = d,
            start = start,
            formula = flist,
            formula_parsed = fp )
    }
    
    return( result )
    
}

# EXAMPLES
if ( FALSE ) {

library(rethinking)

f <- list(
    y ~ dnorm(mu,sigma),
    mu ~ a + aj + (b+bj)*x,
    c(aj,bj)|id ~ dmvnorm( 0 , Sigma_id ),
    a ~ dnorm(0,10),
    b ~ dnorm(0,1),
    sigma ~ dcauchy(0,1)
)

f2 <- list(
    y ~ dnorm(mu,sigma),
    mu ~ a + aj + b*x,
    aj|id ~ dnorm( 0 , sigma_a ),
    a ~ dnorm(0,10),
    b ~ dnorm(0,1),
    sigma ~ dcauchy(0,1),
    sigma_a ~ dcauchy(0,1)
)

f3 <- list(
    y ~ dbinom(1,theta),
    logit(theta) ~ a + aj + b*x,
    aj|id ~ dnorm( 0 , sigma_a ),
    a ~ dnorm(0,10),
    b ~ dnorm(0,1),
    sigma_a ~ dcauchy(0,1)
)

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

f5 <- list(
    y ~ dnorm(mu,sigma),
    mu ~ aj + bj*x,
    c(aj,bj)|id ~ dmvnorm( c(a,b) , Sigma_id ),
    a ~ dnorm(0,10),
    b ~ dnorm(0,1),
    sigma ~ dcauchy(0,1)
)

# simulate data
library(MASS)
N <- 1000 # 1000 cases
J <- 100 # 100 clusters
J2 <- 20
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

m <- map2stan( f2 , data=list(y=y,x=x,id=id) , start=list(a=10,b=0,sigma=3) , sample=TRUE , debug=FALSE )

m2 <- map2stan( f , data=list(y=y,x=x,id=id) , start=list(a=10,b=0,sigma=3) , sample=TRUE , debug=FALSE )


} #EXAMPLES
