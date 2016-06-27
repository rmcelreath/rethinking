# map2stan2
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
        stan_suffix = "_lpdf",
        num_pars = 2,
        par_names = c("mu","sigma"),
        par_bounds = c("","<lower=0>"),
        par_types = c("real","real"),
        out_type = "real",
        par_map = function(k,e,...) {
            # get constraints and add <lower=0> for sigma vector
            constr_list <- get( "constraints" , envir=e )
            sigma_name <- as.character( k[[2]] )
            if ( is.null(constr_list[[sigma_name]]) ) {
                constr_list[[sigma_name]] <- "lower=0"
                assign( "constraints" , constr_list , envir=e )
            }
            return(k);
        },
        vectorized = TRUE
    ),
    LogGaussian = list(
        name = "LogGaussian",
        R_name = "dlnorm",
        stan_name = "lognormal",
        stan_suffix = "_lpdf",
        num_pars = 2,
        par_names = c("mu","sigma"),
        par_bounds = c("","<lower=0>"),
        par_types = c("real","real"),
        out_type = "real",
        par_map = function(k,e,...) {
            # get constraints and add <lower=0> for sigma vector
            constr_list <- get( "constraints" , envir=e )
            sigma_name <- as.character( k[[2]] )
            if ( is.null(constr_list[[sigma_name]]) ) {
                constr_list[[sigma_name]] <- "lower=0"
                assign( "constraints" , constr_list , envir=e )
            }
            return(k);
        },
        vectorized = TRUE
    ),
    StudentT = list(
        name = "StudentT",
        R_name = "dstudent",
        stan_name = "student_t",
        stan_suffix = "_lpdf",
        num_pars = 3,
        par_names = c("nu","mu","sigma"),
        par_bounds = c("<lower=0>","","<lower=0>"),
        par_types = c("real","real","real"),
        out_type = "real",
        par_map = function(k,e,...) {
            return(k);
        },
        vectorized = TRUE
    ),
    MVGaussian = list(
        name = "MVGaussian",
        R_name = "dmvnorm",
        stan_name = "multi_normal",
        stan_suffix = "_lpdf",
        num_pars = 2,
        par_names = c("Mu","Sigma"),
        par_bounds = c("",""),
        par_types = c("vector","cov_matrix"),
        out_type = "vector",
        par_map = function(k,parenvir,n,...) {
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
                    n <- length(k[[1]])-1
                    for ( i in 2:(n+1) ) pars[[i-1]] <- as.character(k[[1]][[i]])
                    vname <- concat( "Mu_" , paste(pars,collapse="") )
                    
                    # get tpars from parent
                    m_tpars1 <- get( "m_tpars1" , envir=parenvir )
                    m_tpars2 <- get( "m_tpars2" , envir=parenvir )
                    
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
        stan_suffix = "_lpdf",
        num_pars = 3,
        par_names = c("Mu","Sigma","Rho"),
        par_bounds = c("","lower=0",""),
        par_types = c("vector","vector","corr_matrix"),
        out_type = "vector",
        par_map = function(k,e,n,...) {
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
            m_tpars2 <- concat( m_tpars2 , indent , Cov_name , " <- quad_form_diag(" , Rho_name , "," , Sigma_name , ");\n" )
            
            # now replace name
            kout[[2]] <- Cov_name
            
            # assign tpars text to parent environment
            assign( "m_tpars1" , m_tpars1 , envir=e )
            assign( "m_tpars2" , m_tpars2 , envir=e )
            
            # get types list and add corr_matrix
            types_list <- get( "types" , envir=e )
            types_list[[Rho_name]] <- concat("corr_matrix")
            assign( "types" , types_list , envir=e )
            
            # get constraints and add <lower=0> for sigma vector
            constr_list <- get( "constraints" , envir=e )
            isnum <- function(x) {
                xn <- suppressWarnings(as.numeric(x))
                return( ifelse(is.na(xn),FALSE,TRUE) ) 
            }
            if ( !isnum(Sigma_name) )
                if ( is.null(constr_list[[Sigma_name]]) ) {
                    constr_list[[Sigma_name]] <- "lower=0"
                    assign( "constraints" , constr_list , envir=e )
                }
            
            # result
            return(kout);
        },
        vectorized = FALSE
    ),
    MVGaussianLNC = list(
        # MVGauss that uses Cholesky factor for correlation matrix
        # also uses non-centered parameterization internally
        # this means substitutes to_vector(z) ~ normal(0,1) for sampling
        # then transforms to correlated samples with:
        #     v <- (diag_pre_multiply(sigma,L_Rho) * z)'
        # but returns usual Rho and sigma parameters for inference
        name = "MVGaussianLNC",
        R_name = "dmvnormNC",
        stan_name = "normal",  # univariate, because uses Cholesky trick
        stan_suffix = "_lpdf",
        num_pars = 2,
        par_names = c("Sigma","Rho"),
        par_bounds = c("lower=0",""),
        par_types = c("vector","corr_matrix"),
        out_type = "vector",
        out_map = function(kout,kin,N,N_txt,e,...) {
            # this function converts c(a,b) style parameters
            # need : to_vector(z) as output
            # but also add transformed parameters to convert back from z
            
            Sigma_name <- as.character(kin[[1]])
            Rho_name <- as.character(kin[[2]])
            L_name <- concat( "L_" , Rho_name )
            
            # build name of z-score matrix
            # dims will be [ num_effects , num_clusters ]
            #z_name <- paste( kout , collapse="" )
            z_name <- concat( "z_" , N_txt )
            
            # add z matrix to parameters
            start <- get( "start" , envir=e )
            start[[z_name]] <- matrix(0,nrow=length(kout),ncol=N)
            assign( "start" , start , envir=e )
            
            # hide z matrix from pars returned in samples
            #pars_hide <- get( "pars_hide" , envir=parent.frame() )
            #pars_hide[[z_name]] <- 1
            #assign( "pars_hide" , pars_hide , envir=parent.frame() )
            
            # coerce type to plain matrix
            types <- get( "types" , envir=parent.frame() )
            types[[z_name]] <- c( "matrix" , concat("[",length(kout),",",N_txt,"]") )
            assign( "types" , types , envir=parent.frame() )
            
            # add v matrix to transformed parameters
            # matrix[N_groups,N_effects] v; 
            # v <- (diag_pre_multiply(sigma,L_Rho) * z)';
            transpars <- get( "transpars" , envir=parent.frame() )
            start <- get( "start" , envir=parent.frame() )
            startp <- get( "start_prior" , envir=parent.frame() )
            v_name <- concat( "v_" , N_txt )
            transpars[[v_name]] <- c(
                concat("matrix[",N_txt,",",length(kout),"] ",v_name),
                concat(v_name," <- (diag_pre_multiply(",Sigma_name,",",L_name,")*",z_name,")'")
            )
            assign( "transpars" , transpars , envir=parent.frame() )
            
            # add conversion v -> c(a,b) in transformed parameters
            # can use transpars list in parent environment to signal
            #   to place parameters in transformed parameters,
            #   with specified code
            pars_elect <- get( "pars_elect" , envir=parent.frame() )
            for ( i in 1:length(kout) ) {
                # clear from start list, so not in parameters block
                apar <- as.character(kout[[i]])
                start[[apar]] <- NULL
                for ( ch in 1:length(startp) ) startp[[ch]][[apar]] <- NULL
                # pattern:
                # apar <- col(v,i) 
                trans_txt <- c(
                    concat( "vector[" , N_txt , "] " , apar ) ,
                    concat( apar , " <- col(" , v_name , "," , i , ")" ) )
                transpars[[apar]] <- trans_txt
                pars_elect[[apar]] <- 1
            }
            assign( "transpars" , transpars , envir=parent.frame() )
            assign( "start" , start , envir=parent.frame() )
            assign( "start_prior" , startp , envir=parent.frame() )
            assign( "pars_elect" , pars_elect , envir=parent.frame() )
            
            # return left hand text for sampling statement
            return( concat("to_vector(",z_name,")") );
        },
        par_map = function(k,e,n_clusters,N_txt,...) {
            # k is list of input parameters
            # n is dimension of multi_normal
            # e is calling environment
            
            kout <- list('0','1') # z ~ normal(0,1)
            indent <- "    "
            
            ###########
            # Sigma and Rho
            
            Sigma_name <- as.character(k[[1]])
            Rho_name <- as.character(k[[2]])
            L_name <- concat( "L_" , Rho_name )
            
            # add Cholesky factor to parameters block
            # use start and types lists
            start <- get( "start" , envir=parent.frame() )
            start[[L_name]] <- diag(n_clusters)
            assign( "start" , start , envir=parent.frame() )
            
            types <- get( "types" , envir=parent.frame() )
            types[[L_name]] <- c("cholesky_factor_corr" , concat("[",n_clusters,"]") )
            assign( "types" , types , envir=parent.frame() )
            
            # convert prior for Rho to Cholesky prior
            # assume using lkj_corr prior, otherwise going to explode
            # also prior for Rho needs to come after this prior,
            #   so is already in m_model_txt at this point,
            #   because parsing is bottom-to-top
            m_model_txt <- get( "m_model_txt" , envir=parent.frame() )
            m_model_txt <- gsub( Rho_name , L_name , m_model_txt , fixed=TRUE )
            m_model_txt <- gsub( "lkj_corr(" , "lkj_corr_cholesky(" , m_model_txt , fixed=TRUE )
            assign( "m_model_txt" , m_model_txt , envir=parent.frame() )
            
            # remove Rho from parameters block 
            start[[Rho_name]] <- NULL
            assign( "start" , start , envir=parent.frame() )
            # have to check start_prior too
            startp <- get( "start_prior" , envir=parent.frame() )
            for ( ch in 1:length(startp) ) startp[[ch]][[Rho_name]] <- NULL
            assign( "start_prior" , startp , envir=parent.frame() )
            
            # need to change this from transpars to GQ
            # intent: add Rho to GENERATED QUANTITITIES
            # do not want it in transformed pars, bc less efficient
            # just need one eval per sample
            # Rho <- L_Rho * L_Rho'
            transpars <- get( "transpars" , envir=parent.frame() )
            transpars[[Rho_name]] <- c(
                concat("matrix[",n_clusters,",",n_clusters,"] ",Rho_name),
                concat(Rho_name," <- ",L_name," * ",L_name,"'")
            )
            assign( "transpars" , transpars , envir=parent.frame() )
            
            # and make sure Rho is returned in samples
            pars_elect <- get( "pars_elect" , envir=parent.frame() )
            pars_elect[[Rho_name]] <- 1
            assign( "pars_elect" , pars_elect , envir=parent.frame() )
            
            # get constraints and add <lower=0> for sigma vector
            constr_list <- get( "constraints" , envir=parent.frame() )
            isnum <- function(x) {
                xn <- suppressWarnings(as.numeric(x))
                return( ifelse(is.na(xn),FALSE,TRUE) ) 
            }
            if ( !isnum(Sigma_name) )
                if ( is.null(constr_list[[Sigma_name]]) ) {
                    constr_list[[Sigma_name]] <- "lower=0"
                    assign( "constraints" , constr_list , envir=e )
                }
            
            # result
            return(kout);
        },
        vectorized = FALSE
    ),
    GaussianProcess = list(
        name = "GaussianProcess",
        R_name = "GPL2",
        stan_name = "multi_normal",
        stan_suffix = "_lpdf",
        num_pars = 4,
        par_names = c("d","eta_sq","rho_sq","sig_sq"),
        par_bounds = c("","<lower=0>","<lower=0>","<lower=0>"),
        par_types = c("matrix","real","real","real"),
        out_type = "vector",
        par_map = function(k,e,n,N_txt,...) {
            indent <- "    "
            kout <- list()
            # need to fill Mu with zeros
            # make sure Mu matches length
            kout[[1]] <- concat( "rep_vector(0," , N_txt , ")" )
            
            # need to construct Sigma
            Cov_name <- concat( "SIGMA_" , k[[1]] )
            kout[[2]] <- Cov_name
            
            ## handle temp cov matrix declaration
            # get model declarations from parent
            m_tmp <- get( "m_model_declare" , envir=e )
            # add local covariance matrix variable
            # will then build it from 'd' matrix in model block
            m_tmp <- concat( m_tmp , indent , "matrix[",N_txt,",",N_txt,"] " , Cov_name , ";\n" )
            # assign declarations text to parent environment
            assign( "m_model_declare" , m_tmp , envir=e )
            
            ## handle loop to construct cov matrix in model block
            # do this in 'm_model_txt' so is just before the ~ statement
            m_tmp <- get( "m_model_txt" , envir=e )
            m_tmp <- concat( m_tmp , indent , "for ( i in 1:(" , N_txt , "-1) )\n" )
            m_tmp <- concat( m_tmp , indent , indent , "for ( j in (i+1):" , N_txt , " ) {\n" )
            m_tmp <- concat( m_tmp , indent , indent , indent , Cov_name , "[i,j] <- " , k[[2]] , "*exp(-" , k[[3]] , "*pow(" , k[[1]] , "[i,j],2));\n" )
            m_tmp <- concat( m_tmp , indent , indent , indent , Cov_name , "[j,i] <- " , Cov_name , "[i,j];\n" )
            m_tmp <- concat( m_tmp , indent , indent , "}\n" )
            m_tmp <- concat( m_tmp , indent , "for ( k in 1:" , N_txt , " )\n" )
            m_tmp <- concat( m_tmp , indent , indent , Cov_name ,"[k,k] <- " , k[[2]] , " + " , k[[4]] , ";\n" )
            assign( "m_model_txt" , m_tmp , envir=e )
            
            # have to make sure the distance matrix is declared in the data block
            # so add it to used_predictors master list
            # this also let's us keep text dimension declaration
            fp_tmp <- get( "fp" , envir=e )
            var_name <- as.character(k[[1]])
            fp_tmp[['used_predictors']][[ var_name ]] <- list( var=var_name , N=c(N_txt,N_txt) , type="matrix" )
            assign( "fp" , fp_tmp , envir=e )
            
            # make sure parameters have positive contraint
            constr_list <- get( "constraints" , envir=e )
            isnum <- function(x) {
                xn <- suppressWarnings(as.numeric(x))
                return( ifelse(is.na(xn),FALSE,TRUE) ) 
            }
            eta_name <- as.character( k[[2]] )
            rho_name <- as.character( k[[3]] )
            sig_name <- as.character( k[[4]] )
            if ( !isnum(eta_name) )
                if ( is.null(constr_list[[ eta_name ]]) ) 
                    constr_list[[eta_name]] <- "lower=0"
            if ( !isnum(rho_name) )
                if ( is.null(constr_list[[ rho_name ]]) ) 
                    constr_list[[rho_name]] <- "lower=0"
            if ( !isnum(sig_name) )
                if ( is.null(constr_list[[ sig_name ]]) ) 
                    constr_list[[sig_name]] <- "lower=0"
            assign( "constraints" , constr_list , envir=e )
            
            # return 2 element parameter vector for multi_normal
            return(kout);
        },
        vectorized = FALSE
    ),
    Cauchy = list(
        name = "Cauchy",
        R_name = "dcauchy",
        stan_name = "cauchy",
        stan_suffix = "_lpdf",
        num_pars = 2,
        par_names = c("location","scale"),
        par_bounds = c("","<lower=0>"),
        par_types = c("real","real"),
        out_type = "real",
        par_map = function(k,e,...) {
            # get constraints and add <lower=0> for scale
            constr_list <- get( "constraints" , envir=e )
            scale_name <- as.character( k[[2]] )
            if ( is.null(constr_list[[scale_name]]) ) {
                constr_list[[scale_name]] <- "lower=0"
                assign( "constraints" , constr_list , envir=e )
            }
            return(k);
        },
        vectorized = TRUE
    ),
    Ordered = list(
        name = "Ordered",
        R_name = "dordlogit",
        stan_name = "ordered_logistic",
        stan_suffix = "_lpmf",
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
                #type_list[[ cuts_name ]] <- concat( "ordered[" , length(cutpars) , "]" )
                type_list[[ cuts_name ]] <- "ordered"
                assign( "types" , type_list , e )
                # rename
                k[[2]] <- cuts_name
            } else {
                # just a name, hopefully
                cuts_name <- as.character(k[[2]])
                # for now, user must specify types=list(cutpoints="ordered")
                # check
                type_list <- get( "types" , e )
                if ( is.null(type_list[[cuts_name]]) ) {
                    #message(concat("Warning: No explicit type declared for ",cuts_name))
                    type_list[[ cuts_name ]] <- "ordered"
                    assign( "types" , type_list , e )
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
        stan_suffix = "_lpdf",
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
        stan_suffix = "_lpdf",
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
        stan_suffix = "_lpdf",
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
        stan_suffix = "_lpdf",
        num_pars = 2,
        par_names = c("min","max"),
        par_bounds = c("",""),
        par_types = c("real","real"),
        out_type = "real",
        out_map = function(kout,kin,...) {
            # need to coerce bounded constraints on the parameter
            constr_list <- get( "constraints" , envir=parent.frame() )
            kout <- as.character(kout)
            if ( is.null(constr_list[[kout]]) ) {
                # added concat() here instead of as.character to handle "-1" parsing as c("-","1")
                # same problem might exist in other templates?
                constr_list[[kout]] <- concat( "lower=" , concat(kin[[1]]) , ",upper=" , concat(kin[[2]]) )
                assign( "constraints" , constr_list , envir=parent.frame() )
            }
            # add comment tag in front, so that uniform dist isn't computed during sims
            # the bounded contraints effectively establish uniform prior in Stan
            kout <- concat( "// " , kout )
            return(kout)
        },
        par_map = function(k,...) {
            return(k);
        },
        vectorized = TRUE
    ),
    Binomial = list(
        name = "Binomial",
        R_name = "dbinom",
        stan_name = "binomial",
        stan_suffix = "_lpmf",
        nat_link = "logit",
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
    Geometric1 = list(
        name = "Geometric1",
        R_name = "dgeom",
        stan_name = "increment_log_prob",
        stan_code = 
"target += (bernoulli_lpmf(1|PAR1) + OUTCOME*bernoulli_lpmf(0|PAR1));",
        stan_dev = 
"dev <- dev + (-2)*(bernoulli_lpmf(1|PAR1) + OUTCOME*bernoulli_lpmf(0|PAR1));",
        num_pars = 3,
        par_names = c("prob"),
        par_bounds = c("<lower=0,upper=1>"),
        par_types = c("real"),
        out_type = "int",
        par_map = function(k,...) {
            return(k);
        },
        vectorized = FALSE
    ),
    Beta = list(
        name = "Beta",
        R_name = "dbeta",
        stan_name = "beta",
        stan_suffix = "_lpdf",
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
    Beta2 = list(
        name = "Beta2",
        R_name = "dbeta2",
        stan_name = "beta",
        stan_suffix = "_lpdf",
        num_pars = 2,
        par_names = c("prob","theta"),
        par_bounds = c("<lower=0,upper=1>","<lower=0>"),
        par_types = c("real","real"),
        out_type = "real",
        par_map = function(k,...) {
            p_name <- k[[1]];
            theta_name <- as.character(k[[2]]);
            k[[1]] <- concat(p_name,"*",theta_name);
            k[[2]] <- concat("(1-",p_name,")*",theta_name);

            # need to coerce bounded constraints on theta parameter
            constr_list <- get( "constraints" , envir=parent.frame() )
            if ( is.null(constr_list[[theta_name]]) ) {
                constr_list[[theta_name]] <- "lower=0"
                assign( "constraints" , constr_list , envir=parent.frame() )
            }

            return(k);
        },
        vectorized = TRUE
    ),
    BetaBinomial = list(
        name = "BetaBinomial",
        R_name = "dbetabinom",
        stan_name = "beta_binomial",
        stan_suffix = "_lpmf",
        num_pars = 3,
        par_names = c("size","prob","theta"),
        par_bounds = c("<lower=1>","<lower=0>","<lower=0>"),
        par_types = c("int","real","real"),
        out_type = "int",
        par_map = function(k,...) {
            p_name <- k[[2]];
            theta_name <- as.character(k[[3]]);
            k[[2]] <- concat(p_name,"*",theta_name);
            k[[3]] <- concat("(1-",p_name,")*",theta_name);

            # need to coerce bounded constraints on theta parameter
            constr_list <- get( "constraints" , envir=parent.frame() )
            if ( is.null(constr_list[[theta_name]]) ) {
                constr_list[[theta_name]] <- "lower=0"
                assign( "constraints" , constr_list , envir=parent.frame() )
            }

            return(k);
        },
        vectorized = TRUE
    ),
    Poisson = list(
        name = "Poisson",
        R_name = "dpois",
        stan_name = "poisson",
        stan_suffix = "_lpmf",
        nat_link = "log",
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
        stan_suffix = "_lpdf",
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
        stan_suffix = "_lpdf",
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
        stan_suffix = "_lpdf",
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

            # need to coerce bounded constraints on scale parameter
            constr_list <- get( "constraints" , envir=parent.frame() )
            if ( is.null(constr_list[[scale]]) ) {
                constr_list[[scale]] <- "lower=0"
                assign( "constraints" , constr_list , envir=parent.frame() )
            }
            
            return(k);
        },
        vectorized = TRUE
    ),
    GammaPoisson = list(
        name = "GammaPoisson",
        R_name = "dgampois",
        stan_name = "neg_binomial_2",
        stan_suffix = "_lpmf",
        num_pars = 2,
        par_names = c("alpha","beta"),
        par_bounds = c("<lower=0>","<lower=0>"),
        par_types = c("real","real"),
        out_type = "int",
        par_map = function(k,...) {
            mu_name <- k[[1]];
            scale_name <- as.character(k[[2]]);
            k[[1]] <- concat(mu_name);
            k[[2]] <- concat(scale_name);

            # need to coerce bounded constraints on scale parameter
            constr_list <- get( "constraints" , envir=parent.frame() )
            if ( is.null(constr_list[[scale_name]]) ) {
                constr_list[[scale_name]] <- "lower=0"
                assign( "constraints" , constr_list , envir=parent.frame() )
            }

            return(k);
        },
        vectorized = TRUE
    ),
    ZeroAugmentedGamma2 = list(
        name = "ZeroAugmentedGamma2",
        R_name = "dzagamma2",
        stan_name = "increment_log_prob",
        stan_code = 
"if (OUTCOME == 0)
target += (bernoulli_lpmf(1|PAR1));
else
target += (bernoulli_lpmf(0|PAR1) + gamma_lpdf(OUTCOME|PAR2,PAR3));",
        stan_dev = 
"if (OUTCOME == 0)
dev <- dev + (-2)*(bernoulli_lpmf(1|PAR1));
else
dev <- dev + (-2)*(bernoulli_lpmf(0|PAR1) + gamma_lpdf(OUTCOME|PAR2,PAR3));",
        num_pars = 3,
        par_names = c("prob","alpha","beta"),
        par_bounds = c("<lower=0,upper=1>","<lower=0>","lower=0"),
        par_types = c("real","real","real"),
        out_type = "real",
        par_map = function(k,...) {
            # alpha is shape
            # beta is rate
            # mu = alpha/scale
            # scale = 1/beta
            mu <- as.character(k[[2]])
            scale <- as.character(k[[3]])
            k[[2]] <- concat( mu , "[i]/" , scale )
            k[[3]] <- concat( "1/" , scale )

            # need to coerce bounded constraints on scale parameter
            constr_list <- get( "constraints" , envir=parent.frame() )
            if ( is.null(constr_list[[scale]]) ) {
                constr_list[[scale]] <- "lower=0"
                assign( "constraints" , constr_list , envir=parent.frame() )
            }

            return(k);
        },
        vectorized = FALSE
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
      target += (log_sum_exp(bernoulli_lpmf(1|PAR1),
        bernoulli_lpmf(0|PAR1) + poisson_lpmf(OUTCOME|PAR2)));
    else
      target += (bernoulli_lpmf(0|PAR1) + poisson_lpmf(OUTCOME|PAR2));",
        stan_dev = 
"if (OUTCOME == 0)
      dev <- dev + (-2)*(log_sum_exp(bernoulli_lpmf(1|PAR1),
        bernoulli_lpmf(0|PAR1) + poisson_lpmf(OUTCOME|PAR2)));
    else
      dev <- dev + (-2)*(bernoulli_lpmf(0|PAR1) + poisson_lpmf(OUTCOME|PAR2));",
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
target += (log_sum_exp(bernoulli_lpmf(1|PAR1),
    bernoulli_lpmf(0|PAR1) + binomial_lpmf(OUTCOME|PAR2,PAR3)));
else
target += (bernoulli_lpmf(0|PAR1) + binomial_lpmf(OUTCOME|PAR2,PAR3));",
        stan_dev = 
"if (OUTCOME == 0)
dev <- dev + (-2)*(log_sum_exp(bernoulli_lpmf(1|PAR1),
    bernoulli_lpmf(0|PAR1) + binomial_lpmf(OUTCOME|PAR2,PAR3)));
else
dev <- dev + (-2)*(bernoulli_lpmf(0|PAR1) + binomial_lpmf(OUTCOME|PAR2,PAR3));",
        num_pars = 2,
        par_names = c("prob1","size","prob2"),
        par_bounds = c("",""),
        par_types = c("real","int","real"),
        out_type = "int",
        par_map = function(k,...) {
            return(k);
        },
        vectorized = FALSE
    ),
    Categorical1 = list(
        name = "Categorical1",
        R_name = "dcategorical",
        stan_name = "increment_log_prob",
        stan_code = 
"{
        vector[PAR1] theta;
        PAR2 
        OUTCOME ~ categorical( softmax(theta) );
    }",
        stan_dev = 
"{
        vector[PAR1] theta;
        PAR2 
        dev <- dev + (-2)*categorical_lpmf( OUTCOME | softmax(theta) );
    }",
        num_pars = 1,
        par_names = c("prob"),
        par_bounds = c("<lower=0,upper=1>"),
        par_types = c("real"),
        out_type = "int",
        par_map = function(k,...) {
            # k should be softmax call
            # convert to list of names with indices
            new_k <- as.character(k[[1]])
            # add length to front, overwriting 'softmax' function name
            num_scores <- length(new_k)-1
            new_k[1] <- num_scores
            # now need to build the code that populates the theta vector of inputs to softmax
            theta_txt <- ""
            for ( i in 1:num_scores ) {
                PAR <- new_k[i+1]
                if ( is.na(as.numeric(PAR)) ) PAR <- concat(PAR,"[i]")
                theta_txt <- concat( theta_txt , "theta[" , i , "] <- " , PAR , ";" )
                if ( i < num_scores ) theta_txt <- concat( theta_txt , "\n        " )
            }
            return(c(new_k[1],theta_txt));
        },
        vectorized = FALSE
    )
)

