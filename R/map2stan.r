# map2stan2 - rewrite of compilation algorithm
# use templates this time
# build all parts of stan code at same time, as pass through formulas

# template design:
# templates map R density functions onto Stan functions

# to-do:
# (*) need to do p*theta and (1-p)*theta multiplication outside likelihood in betabinomial (similar story for gammapoisson) --- or is there an operator for pairwise vector multiplication?
# (*) nobs calculation after fitting needs to account for aggregated binomial
# (*) fix detection of correlation matrix type, "rho" versus "Rho" -- should just insert into type list at dmvnorm2
# (-) handle improper input more gracefully
# (-) add "as is" formula type, with quoted text on RHS to dump into Stan code

########################################
# distribution function templates

concat <- function( ... ) {
    paste( ... , collapse="" , sep="" )
}

map2stan.templates <- list(
    Gaussian = list(
        name = "Gaussian",
        R_name = "dnorm",
        stan_name = "normal",
        num_pars = 2,
        par_names = c("mu","sigma"),
        par_bounds = c("","<lower=0>"),
        par_types = c("real","real"),
        out_type = "real",
        par_map = function(k,...) {
            return(k);
        },
        vectorized = TRUE
    ),
    StudentT = list(
        name = "StudentT",
        R_name = "dstudent",
        stan_name = "student_t",
        num_pars = 3,
        par_names = c("nu","mu","sigma"),
        par_bounds = c("<lower=0>","","<lower=0>"),
        par_types = c("real","real","real"),
        out_type = "real",
        par_map = function(k,...) {
            return(k);
        },
        vectorized = TRUE
    ),
    MVGaussian = list(
        name = "MVGaussian",
        R_name = "dmvnorm",
        stan_name = "multi_normal",
        num_pars = 2,
        par_names = c("Mu","Sigma"),
        par_bounds = c("",""),
        par_types = c("vector","cov_matrix"),
        out_type = "vector",
        par_map = function(k,n,parenvir,...) {
            # k is list of input parameters
            # n is dimension of multi_normal
            kout <- k
            
            # make sure Mu matches length
            if ( class(k[[1]])=="numeric" ) {
                # numeric means -- make sure match dimension
                kout[[1]] <- concat( "rep_vector(" , k[[1]] , "," , n , ")" )
            }
            
            indent <- "    "
            
            # check for vector of parameters in Mu
            if ( class(k[[1]])=="call" ) {
                fname <- as.character(k[[1]][[1]])
                if ( fname=="c" ) {
                    # vector, so constuct vector data type in Stan code
                    # and add vector to transformed parameters
                    pars <- list()
                    for ( i in 2:(n+1) ) pars[[i-1]] <- as.character(k[[1]][[i]])
                    vname <- concat( "Mu_" , paste(pars,collapse="") )
                    
                    # add declaration to transformed parameters
                    m_tpars1 <- concat( m_tpars1 , indent , "vector[" , n , "] " , vname , ";\n" )
                    # add transformation
                    m_tpars2 <- concat( m_tpars2 , indent , "for ( j in 1:" , n , " ) {\n" )
                    for ( i in 1:n ) {
                        m_tpars2 <- concat( m_tpars2 , indent,indent , vname , "[" , i , "] <- " , pars[[i]] , ";\n" )
                    }
                    m_tpars2 <- concat( m_tpars2 , indent , "}\n" )
                    
                    # assign tpars text to parent environment?
                    assign( "m_tpars1" , m_tpars1 , envir=parenvir )
                    assign( "m_tpars2" , m_tpars2 , envir=parenvir )
                    
                    # finally, replace k[[1]] with vector name
                    kout[[1]] <- vname
                }
            }
            
            # result
            return(kout);
        },
        vectorized = FALSE
    ),
    MVGaussianSRS = list(
        # MVGauss with diag(sigma)*Rho*diag(sigma) for covariance
        name = "MVGaussianSRS",
        R_name = "dmvnorm2",
        stan_name = "multi_normal",
        num_pars = 3,
        par_names = c("Mu","Sigma","Rho"),
        par_bounds = c("","<lower=0>",""),
        par_types = c("vector","vector","cor_matrix"),
        out_type = "vector",
        par_map = function(k,n,e,...) {
            # k is list of input parameters
            # n is dimension of multi_normal
            # e is calling environment
            
            # only going to need two slots in result
            kout <- k
            kout[[3]] <- NULL
            
            ###########
            # Mu
            
            indent <- "    "
            
            # make sure Mu matches length
            if ( class(k[[1]])=="numeric" ) {
                # numeric means -- make sure match dimension
                kout[[1]] <- concat( "rep_vector(" , k[[1]] , "," , n , ")" )
            }
            # check for vector of parameters in Mu
            if ( class(k[[1]])=="call" ) {
                fname <- as.character(k[[1]][[1]])
                if ( fname=="c" ) {
                    # vector, so constuct vector data type in Stan code
                    # and add vector to transformed parameters
                    pars <- list()
                    for ( i in 2:(n+1) ) pars[[i-1]] <- as.character(k[[1]][[i]])
                    vname <- concat( "Mu_" , paste(pars,collapse="") )
                    
                    # get tpars from parent
                    m_tpars1 <- get( "m_tpars1" , envir=e )
                    m_tpars2 <- get( "m_tpars2" , envir=e )
                    
                    # add declaration to transformed parameters
                    m_tpars1 <- concat( m_tpars1 , indent , "vector[" , n , "] " , vname , ";\n" )
                    # add transformation
                    m_tpars2 <- concat( m_tpars2 , indent , "for ( j in 1:" , n , " ) {\n" )
                    for ( i in 1:n ) {
                        m_tpars2 <- concat( m_tpars2 , indent,indent , vname , "[" , i , "] <- " , pars[[i]] , ";\n" )
                    }
                    m_tpars2 <- concat( m_tpars2 , indent , "}\n" )
                    
                    # assign tpars text to parent environment
                    assign( "m_tpars1" , m_tpars1 , envir=e )
                    assign( "m_tpars2" , m_tpars2 , envir=e )
                    
                    # finally, replace k[[1]] with vector name
                    kout[[1]] <- vname
                }
            }
            
            ###########
            # Sigma and Rho
            # construct covariance matrix from
            #   diag_matrix(Sigma)*Rho*diag_matrix(Sigma)
            # Then can specify separate priors on Sigma and Rho
            # need to use transformed parameter for construction, 
            #   so calc not repeated in loop
            # Naming convention for cov_matrix: SRS_SigmaRho
            
            Sigma_name <- as.character(k[[2]])
            Rho_name <- as.character(k[[3]])
            Cov_name <- concat( "SRS_" , Sigma_name , Rho_name )
            
            # get tpars from parent
            m_tpars1 <- get( "m_tpars1" , envir=e )
            m_tpars2 <- get( "m_tpars2" , envir=e )
            
            # build transformed parameter
            m_tpars1 <- concat( m_tpars1 , indent , "cov_matrix[" , n , "] " , Cov_name , ";\n" )
            m_tpars2 <- concat( m_tpars2 , indent , Cov_name , " <- diag_matrix(" , Sigma_name , ")*" , Rho_name , "*diag_matrix(" , Sigma_name , ");\n" )
            
            # now replace name
            kout[[2]] <- Cov_name
            
            # assign tpars text to parent environment
            assign( "m_tpars1" , m_tpars1 , envir=e )
            assign( "m_tpars2" , m_tpars2 , envir=e )
            
            # result
            return(kout);
        },
        vectorized = FALSE
    ),
    Cauchy = list(
        name = "Cauchy",
        R_name = "dcauchy",
        stan_name = "cauchy",
        num_pars = 2,
        par_names = c("location","scale"),
        par_bounds = c("","<lower=0>"),
        par_types = c("real","real"),
        out_type = "real",
        par_map = function(k,...) {
            return(k);
        },
        vectorized = TRUE
    ),
    Ordered = list(
        name = "Ordered",
        R_name = "dordlogit",
        stan_name = "ordered_logistic",
        num_pars = 2,
        par_names = c("eta","cutpoints"),
        par_bounds = c("",""),
        par_types = c("real","vector"),
        out_type = "int",
        par_map = function(k,e,...) {
            # linear model eta doesn't need any special treatment
            # the cutpoints need to be declared as type ordered with length K-1,
            #   where K is number of levels in outcome
            # try to do as much here as we can
            
            # make sure cutpoints are not in R vector form
            cuts <- k[[2]]
            if ( class(cuts)=="call" ) {
                # convert to single name
                cuts_name <- "cutpoints"
                # get start values for individual pars in start list
                s <- get( "start" , e )
                cutpars <- list()
                cutstart <- list()
                for ( i in 2:length(cuts) ) {
                    cutpars[[i-1]] <- as.character(cuts[[i]])
                    cutstart[[i-1]] <- s[[ cutpars[[i-1]] ]]
                    s[[ cutpars[[i-1]] ]] <- NULL # remove old naming
                }
                s[[ cuts_name ]] <- as.numeric(cutstart)
                assign( "start" , s , e )
                # insert explicit type into types list
                type_list <- get( "types" , e )
                type_list[[ cuts_name ]] <- concat( "ordered[" , length(cutpars) , "]" )
                assign( "types" , type_list , e )
                # rename
                k[[2]] <- cuts_name
            } else {
                # just a name, hopefully
                cuts_name <- as.character(k[[2]])
                # for now, user must specify types=list(cutpoints="ordered[K]")
                # check
                type_list <- get( "types" , e )
                if ( is.null(type_list[[cuts_name]]) ) {
                    message(concat("Warning: No explicit type declared for ",cuts_name))
                }
                k[[2]] <- cuts_name
            }
            
            # add [i] to eta -- ordered not vectorized
            # this should work for both linear models and variables
            eta <- as.character(k[[1]])
            k[[1]] <- concat( eta , "[i]" )
            
            # result
            return(k);
        },
        vectorized = FALSE
    ),
    LKJ_Corr = list(
        # LKJ corr_matrix; eta=1 is uniform correlation matrices
        name = "LKJ_Corr",
        R_name = "dlkjcorr",
        stan_name = "lkj_corr",
        num_pars = 1,
        par_names = c("eta"),
        par_bounds = c("<lower=0>"),
        par_types = c("real"),
        out_type = "corr_matrix",
        par_map = function(k,...) {
            return(k);
        },
        vectorized = FALSE
    ),
    invWishart = list(
        # ye olde inverse wishart prior -- nu is df
        name = "invWishart",
        R_name = "dinvwishart",
        stan_name = "inv_wishart",
        num_pars = 2,
        par_names = c("nu","Sigma"),
        par_bounds = c("<lower=0>",""),
        par_types = c("real","matrix"),
        out_type = "matrix",
        par_map = function(k,...) {
            # process Sigma
            # check if diag(n) format and translate to Stan broadcast code
            if ( class(k[[2]])=="call" ) {
                # could be function call
                fname <- as.character(k[[2]][[1]])
                if ( fname=="diag" ) {
                    n <- as.integer(k[[2]][[2]])
                    k[[2]] <- concat("diag_matrix(rep_vector(1,",n,"))")
                }
            }
            # result
            return(k);
        },
        vectorized = FALSE
    ),
    Laplace = list(
        name = "Laplace",
        R_name = "dlaplace",
        stan_name = "double_exponential",
        num_pars = 2,
        par_names = c("location","scale"),
        par_bounds = c("<lower=0>","<lower=0>"),
        par_types = c("real","real"),
        out_type = "real",
        par_map = function(k,...) {
            # Stan uses 1/lambda for second parameter
            lambda <- as.character(k[[2]])
            k[[2]] <- concat( "1/",lambda )
            return(k);
        },
        vectorized = TRUE
    ),
    Uniform = list(
        name = "Uniform",
        R_name = "dunif",
        stan_name = "uniform",
        num_pars = 2,
        par_names = c("alpha","beta"),
        par_bounds = c("<lower=0>","<lower=0>"),
        par_types = c("real","real"),
        out_type = "real",
        par_map = function(k,...) {
            return(k);
        },
        vectorized = TRUE
    ),
    Binomial = list(
        name = "Binomial",
        R_name = "dbinom",
        stan_name = "binomial",
        num_pars = 2,
        par_names = c("size","prob"),
        par_bounds = c("<lower=1>","<lower=0,upper=1>"),
        par_types = c("int","real"),
        out_type = "int",
        par_map = function(k,...) {
            return(k);
        },
        vectorized = TRUE
    ),
    Beta = list(
        name = "Beta",
        R_name = "dbeta",
        stan_name = "beta",
        num_pars = 2,
        par_names = c("alpha","beta"),
        par_bounds = c("<lower=0>","<lower=0>"),
        par_types = c("real","real"),
        out_type = "real",
        par_map = function(k,...) {
            return(k);
        },
        vectorized = TRUE
    ),
    BetaBinomial = list(
        name = "BetaBinomial",
        R_name = "dbetabinom",
        stan_name = "beta_binomial",
        num_pars = 3,
        par_names = c("size","prob","theta"),
        par_bounds = c("<lower=1>","<lower=0>","<lower=0>"),
        par_types = c("int","real","real"),
        out_type = "int",
        par_map = function(k,...) {
            p_name <- k[[2]];
            theta_name <- k[[3]];
            k[[2]] <- concat(p_name,"*",theta_name);
            k[[3]] <- concat("(1-",p_name,")*",theta_name);
            return(k);
        },
        vectorized = TRUE
    ),
    Poisson = list(
        name = "Poisson",
        R_name = "dpois",
        stan_name = "poisson",
        num_pars = 1,
        par_names = c("lambda"),
        par_bounds = c("<lower=0>"),
        par_types = c("real"),
        out_type = "int",
        par_map = function(k,...) {
            return(k);
        },
        vectorized = TRUE
    ),
    Exponential = list(
        name = "Exponential",
        R_name = "dexp",
        stan_name = "exponential",
        num_pars = 1,
        par_names = c("lambda"),
        par_bounds = c("<lower=0>"),
        par_types = c("real"),
        out_type = "real",
        par_map = function(k,...) {
            return(k);
        },
        vectorized = TRUE
    ),
    Gamma = list(
        name = "Gamma",
        R_name = "dgamma",
        stan_name = "gamma",
        num_pars = 2,
        par_names = c("alpha","beta"),
        par_bounds = c("<lower=0>","lower=0"),
        par_types = c("real","real"),
        out_type = "real",
        par_map = function(k,...) {
            # alpha is shape
            # beta is rate
            return(k);
        },
        vectorized = TRUE
    ),
    Gamma2 = list(
        name = "Gamma2",
        R_name = "dgamma2",
        stan_name = "gamma",
        num_pars = 2,
        par_names = c("alpha","beta"),
        par_bounds = c("<lower=0>","lower=0"),
        par_types = c("real","real"),
        out_type = "real",
        par_map = function(k,...) {
            # alpha is shape
            # beta is rate
            # mu = alpha/scale
            # scale = 1/beta
            mu <- as.character(k[[1]])
            scale <- as.character(k[[2]])
            k[[1]] <- concat( mu , "/" , scale )
            k[[2]] <- concat( "1/" , scale )
            return(k);
        },
        vectorized = TRUE
    ),
    GammaPoisson = list(
        name = "GammaPoisson",
        R_name = "dgampois",
        stan_name = "neg_binomial",
        num_pars = 2,
        par_names = c("alpha","beta"),
        par_bounds = c("<lower=0>","<lower=0>"),
        par_types = c("real","real"),
        out_type = "int",
        par_map = function(k,...) {
            mu_name <- k[[1]];
            scale_name <- k[[2]];
            k[[1]] <- concat(mu_name,"/",scale_name);
            k[[2]] <- concat("1/",scale_name);
            return(k);
        },
        vectorized = TRUE
    ),
    ZeroInflatedPoisson = list(
        # not built into Stan, but can build from log_prob calculations
        # need to flag by using 'increment_log_prob' as distribution name
        # then provide separate model{} and gq{} code segments
        name = "ZeroInflatedPoisson",
        R_name = "dzipois",
        stan_name = "increment_log_prob",
        stan_code = 
"if (OUTCOME == 0)
increment_log_prob(log_sum_exp(bernoulli_log(1,PAR1),
    bernoulli_log(0,PAR1) + poisson_log(OUTCOME,PAR2)));
else
increment_log_prob(bernoulli_log(0,PAR1) + poisson_log(OUTCOME,PAR2));",
        stan_dev = 
"if (OUTCOME == 0)
dev <- dev + (-2)*(log_sum_exp(bernoulli_log(1,PAR1),
    bernoulli_log(0,PAR1) + poisson_log(OUTCOME,PAR2)));
else
dev <- dev + (-2)*(bernoulli_log(0,PAR1) + poisson_log(OUTCOME,PAR2));",
        num_pars = 2,
        par_names = c("p","lambda"),
        par_bounds = c("",""),
        par_types = c("real","real"),
        out_type = "int",
        par_map = function(k,...) {
            return(k);
        },
        vectorized = FALSE
    ),
    ZeroInflatedBinomial = list(
        # not built into Stan, but can build from log_prob calculations
        # need to flag by using 'increment_log_prob' as distribution name
        # then provide separate model{} and gq{} code segments
        name = "ZeroInflatedBinomial",
        R_name = "dzibinom",
        stan_name = "increment_log_prob",
        stan_code = 
"if (OUTCOME == 0)
increment_log_prob(log_sum_exp(bernoulli_log(1,PAR1),
    bernoulli_log(0,PAR1) + binomial_log(OUTCOME,PAR2,PAR3)));
else
increment_log_prob(bernoulli_log(0,PAR1) + binomial_log(OUTCOME,PAR2,PAR3));",
        stan_dev = 
"if (OUTCOME == 0)
dev <- dev + (-2)*(log_sum_exp(bernoulli_log(1,PAR1),
    bernoulli_log(0,PAR1) + binomial_log(OUTCOME,PAR2,PAR3)));
else
dev <- dev + (-2)*(bernoulli_log(0,PAR1) + binomial_log(OUTCOME,PAR2,PAR3));",
        num_pars = 2,
        par_names = c("prob1","size","prob2"),
        par_bounds = c("",""),
        par_types = c("real","int","real"),
        out_type = "int",
        par_map = function(k,...) {
            return(k);
        },
        vectorized = FALSE
    )
)


##################
# map2stan itself

map2stan <- function( flist , data , start , pars , constraints=list() , types=list() , sample=TRUE , iter=2000 , chains=1 , debug=FALSE , WAIC=FALSE , ... ) {
    
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
    m_gq <- "" # generated quantities, can build mostly from m_model pieces
    
    ########################################
    # inverse link list
    inverse_links <- list(
        log = 'exp',
        logit = 'inv_logit',
        logodds = 'inv_logit'
    )
    
    templates <- map2stan.templates
    
    ########################################
    # check arguments
    if ( missing(flist) ) stop( "Formula required." )
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
    start.orig <- start
    
    ########################################
    # private functions
    
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
    
    wildpatt <- "[()=*+/ ]"
    
    mygrep <- function( target , replacement , x , add.par=TRUE ) {
        wild <- wildpatt
        pattern <- paste( wild , target , wild , sep="" , collapse="" )
        
        m <- gregexpr( pattern , x )
        
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
    
    # for converting characters not allows by Stan
    undot <- function( astring ) {
        astring <- gsub( "." , "_" , astring , fixed=TRUE )
        astring
    }
    
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
        wild <- wildpatt
        x2 <- paste( " " , x , " " , collapse="" , sep="" )
        pattern <- paste( wild , target , wild , sep="" , collapse="" )
        m <- regexpr( pattern , x2 )
        return(m)
    }
    
    ########################################
    # parse formulas
        
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
            out_type = template$out_type
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
        # find likelihood that contains this lm
        N_name <- "N" # default
        if ( length(fp[['likelihood']]) > 0 ) {
            for ( i in 1:length(fp[['likelihood']]) ) {
                pars <- fp[['likelihood']][[i]][['pars']]
                # find par in likelihood that matches lhs of this lm
                for ( j in 1:length(pars) ) {
                    if ( parameter == pars[[j]] ) {
                        N_name <- fp[['likelihood']][[i]][['N_name']]
                    }
                }
            }
        }
        # result
        list(
            parameter = parameter ,
            RHS = RHS ,
            link = link ,
            N_name = N_name
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
        apar <- as.character( fl[[2]] )
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
    d <- as.list(data)
    
    ##################################
    # check for previous fit object
    # if found, build again to get data right, but use compiled model later
    
    flag_refit <- FALSE
    if ( class(flist)=="map2stan" ) {
        oldfit <- flist
        flist <- oldfit@formula
        flag_refit <- TRUE
    }
    
    ####
    # pass over formulas and extract into correct slots
    # we go from top to bottom, so can see where linear models plug in, when we find them 
    
    for ( i in 1:length(flist) ) {
    
        # test for likelihood
        # if has distribution function on RHS and LHS is variable, then likelihood
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
                    if ( as.character(LHS) %in% names(d) ) {
                        is_likelihood <- TRUE
                    }
                }
            }
        }
        if ( is_likelihood==TRUE ) {
            lik <- extract_likelihood( flist[[i]] )
            fp[['likelihood']] <- listappend( fp[['likelihood']] , lik )
            
            # add outcome to used variables
            fp[['used_predictors']] <- listappend( fp[['used_predictors']] , list( var=undot(lik$outcome) , N=lik$N_name , type=lik$out_type ) )
            
            # check for binomial size variable and mark used
            if ( lik$likelihood=='binomial' | lik$likelihood=='beta_binomial' ) {
                sizename <- as.character(lik$pars[[1]])
                if ( !is.null( d[[sizename]] ) ) {
                    fp[['used_predictors']] <- listappend( fp[['used_predictors']] , list( var=undot(sizename) , N=lik$N_name , type="int" ) )
                }
            }
            if ( lik$template=="ZeroInflatedBinomial" ) {
                # second paramter of zibinom is binomial size variable
                sizename <- as.character(lik$pars[[2]])
                if ( !is.null( d[[sizename]] ) ) {
                    fp[['used_predictors']] <- listappend( fp[['used_predictors']] , list( var=undot(sizename) , N=lik$N_name , type="int" ) )
                }
            }
            next
        }
        
        # test for linear model
        is_linearmodel <- FALSE
        if ( length(flist[[i]][[3]]) > 1 ) {
            fname <- as.character( flist[[i]][[3]][[1]] )
            if ( fname=="+" | fname=="*" | fname=="-" | fname=="/" | fname=="%*%" | fname %in% c('exp','inv_logit') ) {
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
        # can use `|` or `[` to specify varying effects (vector of parameters)
        if ( length( flist[[i]][[2]] ) == 3 ) {
            fname <- as.character( flist[[i]][[2]][[1]] )
            if ( fname=="|" | fname=="[" ) {
                n <- length( fp[['vprior']] )
                fp[['vprior']][[n+1]] <- extract_vprior( flist[[i]] )
                next
            }
        }
        
        # an ordinary prior?
        # check for vectorized LHS
        if ( length( flist[[i]][[2]] ) > 1 ) {
            fname <- as.character( flist[[i]][[2]][[1]] )
            if ( fname=="c" ) {
                # get list of parameters and add each as parsed prior
                np <- length( flist[[i]][[2]] )
                for ( j in 2:np ) {
                    fcopy <- flist[[i]]
                    fcopy[[2]] <- flist[[i]][[2]][[j]]
                    n <- length( fp[['prior']] )
                    fp[['prior']][[n+1]] <- extract_prior( fcopy )
                }
            }
        } else {
            # ordinary simple prior
            n <- length( fp[['prior']] )
            fp[['prior']][[n+1]] <- extract_prior( flist[[i]] )
        }
    }
    
    #####
    # add index brackets in linear models
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
                if ( used > -1 & v %in% names(d) ) {
                    # if variable (not lm), add to used predictors list
                    # nup <- length(fp[['used_predictors']])
                    fp[['used_predictors']][[undot(v)]] <- list( var=undot(v) , N=fp[['lm']][[i]][['N_name']] )
                }
                # add index and undot the name
                fp[['lm']][[i]][['RHS']] <- indicize( v , index , fp[['lm']][[i]][['RHS']] , replace=undot(v) )
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
    
    # pass back through parsed formulas and build Stan code
    
    # priors
    n <- length( fp[['prior']] )
    if ( n > 0 ) {
        for ( i in 1:n ) {
            prior <- fp[['prior']][[i]]
            tmplt <- templates[[prior$template]]
            parstxt <- paste( tmplt$par_map(prior$pars_in) , collapse=" , " )
            txt <- concat( indent , prior$par_out , " ~ " , prior$density , "( " , parstxt , " );" )
            m_model_priors <- concat( m_model_priors , txt , "\n" )
        }#i
    }
    
    # vpriors
    n <- length( fp[['vprior']] )
    if ( n > 0 ) {
        for ( i in 1:n ) {
            vprior <- fp[['vprior']][[i]]
            tmplt <- templates[[vprior$template]]
            N_txt <- concat( "N_" , vprior$group )
            npars <- length(vprior$pars_out)
            
            # lhs -- if vector, need transformed parameter of vector type
            lhstxt <- ""
            if ( length(vprior$pars_out) > 1 ) {
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
            } else {
                # single parameter
                lhstxt <- vprior$pars_out[[1]]
            }
            
            # format parmater inputs
            # use par_map function in template, so ordering etc can change
            rhstxt <- paste( tmplt$par_map( vprior$pars_in , npars , environment() ) , collapse=" , " )
            
            # add text to model code
            m_model_priors <- concat( m_model_priors , indent , "for ( j in 1:" , N_txt , " ) " , lhstxt , "[j] ~ " , vprior$density , "( " , rhstxt , " );\n" )
            
            # declare each parameter with correct type from template
            outtype <- "vector"
            for ( j in 1:length(vprior$pars_out) ) {
                #m_pars <- concat( m_pars , indent , outtype , "[" , N_txt , "] " , vprior$pars_out[[j]] , ";\n" )
            }
            
            # add data declaration for grouping variable
            m_data <- concat( m_data , indent , "int<lower=1> " , N_txt , ";\n" )
            
            # add N count to data
            N <- length( unique( d[[ vprior$group ]] ) )
            d[[ N_txt ]] <- N
            
            # mark grouping variable used
            #fp[['used_predictors']] <- listappend( fp[['used_predictors']] , list(var=vprior$group,N=length(d[[vprior$group]]) ) )
            # check likelihoods for matching length and use that N_name
            if ( length(fp[['likelihood']])>0 ) {
                groupN <- length(d[[vprior$group]])
                for ( j in 1:length(fp[['likelihood']]) ) {
                    if ( fp[['likelihood']][[j]]$N_cases == groupN ) {
                        fp[['used_predictors']] <- listappend( fp[['used_predictors']] , list(var=vprior$group,N=fp[['likelihood']][[j]]$N_name ) )
                    }
                }#j
            } else {
                # just add raw integer length
                fp[['used_predictors']] <- listappend( fp[['used_predictors']] , list(var=vprior$group,N=length(d[[vprior$group]]) ) )
            }
        }#i
    }
    
    if ( debug==TRUE ) { print(m_tpars1); print(m_tpars2); }
    
    # linear models
    # need to do these in reverse (n:1), so intermediate lm's can cascade up
    n <- length( fp[['lm']] )
    if ( n > 0 ) {
        for ( i in n:1 ) {
            linmod <- fp[['lm']][[i]]
            N_txt <- linmod$N_name
            m_model_lm <- concat( m_model_lm , indent , "for ( i in 1:" , N_txt , " ) {\n" )
            m_model_lm <- concat( m_model_lm , indent,indent , linmod$parameter , "[i] <- " , linmod$RHS , ";\n" )
            if ( linmod$link != "identity" ) {
                m_model_lm <- concat( m_model_lm , indent,indent , linmod$parameter , "[i] <- " , inverse_links[[linmod$link]] , "(" , linmod$parameter , "[i]);\n" )
            }
            m_model_lm <- concat( m_model_lm , indent , "}\n" )
            
            # add declaration of linear model local variable
            m_model_declare <- concat( m_model_declare , indent , "vector[" , linmod$N_name , "] " , linmod$parameter , ";\n" )
        }#i
    }
    
    # likelihoods and gq
    
    # build generated quantities
    m_gq <- concat( m_gq , m_model_declare )
    m_gq <- concat( m_gq , indent , "real dev;\n" , indent , "dev <- 0;\n" )
    m_gq <- concat( m_gq , m_model_lm )
    
    # build likelihoods
    n <- length( fp[['likelihood']] )
    if ( n > 0 ) {
        for ( i in 1:n ) {
            lik <- fp[['likelihood']][[i]]
            tmplt <- templates[[lik$template]]
            parstxt_L <- tmplt$par_map( lik$pars , environment() )
            
            # add sampling statement to model block
            outcome <- lik$outcome
            
            if ( tmplt$vectorized==FALSE ) {
                # add loop for non-vectorized distribution
                m_model_lik <- concat( m_model_lik , indent , "for ( i in 1:" , lik$N_name , " )\n" )
                m_gq <- concat( m_gq , indent , "for ( i in 1:" , lik$N_name , " )\n" )
                
                # add [i] to outcome
                outcome <- concat( outcome , "[i]" )
                
                # check for linear model names as parameters and add [i] to each
                if ( length(fp[['lm']])>0 ) {
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
                m_model_lik <- concat( m_model_lik , indent , code_model , "\n" )
                m_gq <- concat( m_gq , indent , code_gq , "\n" )
            } else {
                # regular distribution with ~
                m_model_lik <- concat( m_model_lik , indent , outcome , " ~ " , lik$likelihood , "( " , parstxt , " );\n" )
                m_gq <- concat( m_gq , indent , "dev <- dev + (-2)*" , lik$likelihood , "_log( " , outcome , " , " , parstxt , " );\n" )
            }
            
            # add N variable to data block
            m_data <- concat( m_data , indent , "int<lower=1> " , lik$N_name , ";\n" )
            
            # add number of cases to data list
            d[[ lik$N_name ]] <- as.integer(lik$N_cases)
        }#i
    }
    
    # data from used_predictors list
    n <- length( fp[['used_predictors']] )
    if ( n > 0 ) {
        for ( i in 1:n ) {
            var <- fp[['used_predictors']][[i]]
            type <- "real"
            # integer check
            if ( class(d[[var$var]])=="integer" ) type <- "int"
            # coerce outcome type
            if ( !is.null(var$type) ) type <- var$type
            # build
            m_data <- concat( m_data , indent , type , " " , var$var , "[" , var$N , "];\n" )
        }#i
    }
    
    # declare parameters from start list
    # use any custom constraints in constraints list
    n <- length( start )
    if ( n > 0 ) {
        for ( i in 1:n ) {
            pname <- names(start)[i]
            type <- "real"
            constraint <- ""
            
            if ( class(start[[i]])=="matrix" ) {
                # check for square matrix? just use nrow for now.
                type <- concat( "cov_matrix[" , nrow(start[[i]]) , "]" )
                # corr_matrix check by naming convention
                Rho_check <- grep( "Rho" , pname )
                if ( length(Rho_check)>0 ) type <- concat( "corr_matrix[" , nrow(start[[i]]) , "]" )
            }
            
            # add correct length to numeric vectors (non-matrix)
            if ( type=="real" & length(start[[i]])>1 ) {
                type <- concat( "vector[" , length(start[[i]]) , "]" )
                if ( length(grep("sigma",pname))>0 )
                    type <- concat( "vector<lower=0>[" , length(start[[i]]) , "]" )
                # check for varying effect vector by peeking at vprior list
                # we want symbolic length name here, if varying effect vector
                nvp <- length( fp[['vprior']] )
                if ( nvp > 0 ) {
                    for ( j in 1:nvp ) {
                        pars_out <- fp[['vprior']][[j]]$pars_out
                        for ( k in 1:length(pars_out) ) {
                            if ( pars_out[[k]]==pname ) {
                                type <- concat( "vector[N_" , fp[['vprior']][[j]]$group , "]" )
                            }
                        }#k
                    }#j
                }
            }
            
            # add non-negative restriction to any parameter with 'sigma' in name
            if ( length(grep("sigma",pname))>0 ) {
                if ( type=="real" ) constraint <- "<lower=0>"
            }
            
            # any custom constraint?
            constrainttxt <- constraints[[pname]]
            if ( !is.null(constrainttxt) ) {
                constraint <- concat( "<" , constrainttxt , ">" )
            }
            
            # any custom type?
            mytype <- types[[pname]]
            if ( !is.null(mytype) ) type <- mytype
            
            # add to parameters block
            m_pars <- concat( m_pars , indent , type , constraint , " " , pname , ";\n" )
        }#i
    }
    
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
    
    model_code <- concat( m_data , m_pars , m_tpars1 , m_tpars2 , "model{\n" ,  m_model_declare , m_model_priors , m_model_lm , m_model_lik , "}\n" , m_gq )
    
    if ( debug==TRUE ) cat(model_code)

##############################
# end of Stan code compilation
##############################
    
    ########################################
    # fit model
    
    # build pars vector
    # use ours, unless user provided one
    if ( missing(pars) ) {
        pars <- names(start)
        pars <- c( pars , "dev" )
    }
    
    if ( sample==TRUE ) {
        require(rstan)
        
        # sample
        modname <- deparse( flist[[1]] )
        initlist <- list()
        for ( achain in 1:chains ) initlist[[achain]] <- start
        
        if ( flag_refit==FALSE ) {
            fit <- stan( model_code=model_code , model_name=modname , data=d , init=initlist , iter=iter , chains=chains , pars=pars , ... )
        } else {
            message(concat("Reusing previously compiled model ",oldfit@stanfit@model_name))
            fit <- stan( fit=oldfit@stanfit , model_name=modname , data=d , init=initlist , iter=iter , chains=chains , pars=pars , ... )
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
            names(coef) <- names(start)
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
        fit2 <- stan( fit=fit , init=list(Epost) , data=d , pars="dev" , chains=1 , iter=1 , refresh=-1 )
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
            pars = pars,
            formula = flist,
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
            waic <- WAIC( result , n=0 ) # n=0 to use all available samples
            attr(result,"WAIC") = waic
        }
        
    } else {
        # just return list
        result <- list( 
            call = match.call(), 
            model = model_code,
            data = d,
            start = start,
            pars = pars,
            formula = flist,
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
