# distribution templates for new map2stan - aka ulam

ulam_dists <- list(
    # DEFAULT template used for standardized build function
    # specific distributions can define alternatives
    DEFAULT = list(
        R_name = "DEFAULT",
        Stan_name = "DEFAULT",
        Stan_suffix = "lpdf",
        build = function( left , right , dist_name="XXXX" , suffix="lpdf" , as_log_lik=FALSE , vectorized=TRUE ) {
            right <- as.character( right )
            out <- ""
            for ( j in 1:length(left) ) {
                sym_suffix <- ""

                # check for any multiple choice flags on left symbol
                do_iff <- FALSE
                if( !is.null( attr( left , "iff" ) ) ) {
                    # need to add index to symbols in condition
                    iff_out <- attr( left , "iff" )
                    symlist <- unique( c( left[j] , names(data) , names(symbols) ) )
                    iff_out <- nest_indexes( iff_out , symlist )
                    do_iff <- TRUE
                    vectorized <- FALSE # need explicit loop for conditionals
                }

                # in case of explicit loop, need max index
                # for matrices & arrays, loop over first index only
                max_index <- dim( data[[ left[j] ]] )[1]
                if ( is.null(max_index) ) max_index <- length( data[[ left[j] ]] )

                # outcome var and distribution name
                if ( as_log_lik==FALSE ) {
                    if ( vectorized==FALSE ) {
                        out <- concat( out , indent , "for ( i in 1:" , max_index , " ) " )
                        if ( do_iff == TRUE )
                            out <- concat( out , "\n" , inden(2) , "if ( " , deparse(iff_out) , " ) " )
                        out <- concat( out , left[j] , "[i] ~ " , dist_name , "( " )
                        sym_suffix <- "[i]"
                    } else {
                        out <- concat( out , "    " , left[j] , " ~ " , dist_name , "( " )
                    }
                } 
                # outcome as log_like calc in gq block
                if ( as_log_lik==TRUE ){
                    # log_lik expressions never vectorized, because vectorized distributions return a sum of log_likelihoods, and we need individual terms
                    out <- concat( out , "    " , "for ( i in 1:" , max_index , " ) " )
                    if ( do_iff == TRUE )
                        out <- concat( out , "\n" , inden(2) , "if ( " , deparse(iff_out) , " ) " )
                    out <- concat( out , "log_lik[i] = " , dist_name , "_" , suffix , "( " , left[j] , "[i] | " )
                    sym_suffix <- "[i]"
                }
                # outcome as part of reduce_sum function (multi-threading)
                reduce_context <- FALSE
                if ( as_log_lik=="reduce" ){
                    # for reduce_sum
                    # could be vectorized or not
                    if ( vectorized==TRUE ) {
                        out <- concat( out , "    " , "return " )
                        if ( do_iff == TRUE )
                            out <- concat( out , "\n" , inden(2) , "if ( " , deparse(iff_out) , " ) " )
                        out <- concat( out , dist_name , "_" , suffix , "( " , left[j] , " | " )
                        sym_suffix <- "[start:end]"
                    } else {
                        # something like ordered_logistic
                        # needs loop to sum up terms, because not vectorized
                        out <- concat( out , "   {real lp=0;\n")
                        out <- concat( out , "    " , "for ( i in 1:size(" , left[j] , ") )\n" )
                        out <- concat( out , inden(2) , "lp += " )
                        out <- concat( out , dist_name , "_" , suffix , "( " , left[j] , "[i] | " )
                        #out <- concat( out , inden(1) , "return lp;" )
                        sym_suffix <- "[start:end][i]"
                    }
                    as_log_lik <- TRUE
                    reduce_context <- TRUE
                }

                # add arguments (data or parameters or locals on right-hand side)
                for ( k in 2:length(right) ) {
                    use_sym_suffix <- sym_suffix
                    # parameters are not usually indexed [i] in non-vectorized forms
                    # but data and locals usually are
                    if ( vectorized==FALSE || as_log_lik==TRUE ) {
                        if ( !is.null( symbols[[ right[k] ]] ) ) {
                            # for constants, parameters, and scalar locals, don't add [i]
                            if ( symbols[[ right[k] ]]$type=="par" ) use_sym_suffix <- ""
                            if ( symbols[[ right[k] ]]$type=="local" ) {
                                if ( reduce_context==TRUE & vectorized==FALSE )
                                    # need to index linear model symbol
                                    use_sym_suffix <- "[i]"
                                if ( !is.null(symbols[[ right[k] ]]$dims) ) {
                                    if ( class(symbols[[ right[k] ]]$dims)=="character" | (reduce_context==TRUE & vectorized==TRUE ) )
                                        use_sym_suffix <- ""
                                    else {
                                        # should be a list
                                        if ( symbols[[ right[k] ]]$dims[[1]]=="matrix" ) {
                                            # matrix local, so no loop? have to guess intent is matrix multiplication on right
                                            use_sym_suffix <- ""
                                        }
                                    }
                                }
                            }
                        } else {
                            # not in symbols list. could be a constant or data that hasn't been parsed into symbols yet.
                            if ( !(right[k] %in% names(data)) )
                                use_sym_suffix <- ""
                        }
                    }
                    # attach argument
                    if ( k==length(right) )
                        out <- concat( out , right[k] , use_sym_suffix , " );\n" )
                    else
                        out <- concat( out , right[k] , use_sym_suffix , " , " )
                }#k
                if ( reduce_context==TRUE & vectorized==FALSE ) {
                    # need to close up lp loop context
                    out <- concat( out , inden(1) , "return lp;}\n" )
                }
            }#j
            return( out )
        }
    ),
    uniform = list(
        R_name = "dunif",
        Stan_name = "uniform",
        Stan_suffix = "lpdf",
        pars = 2,
        dims = c( "real" , "real" , "real" ),
        constraints = c( NA , NA , NA ),
        vectorized = TRUE,
        build = function( left , right , as_log_lik=FALSE ) {

            # need to add lower,upper constraints on left variable

            # right: [[1]] dist, [[2]] lower, [[3]] upper
            low <- deparse( right[[2]] )
            upp <- deparse( right[[3]] )
            constraint_txt <- concat( "lower=",low,",upper=",upp )
            left_name <- as.character( left[1] )

            new_symbols <- symbols
            if ( new_symbols[[left_name]]$dims[[1]]=="ordered" )
                warning(concat(left_name," : Constraints not supported for ordered vectors."))
            new_symbols[[left_name]]$constraint <- constraint_txt
            assign( "symbols" , new_symbols , inherits=TRUE )
            
            # call DEFAULT function for rest
            build_func <- distribution_library[["DEFAULT"]]$build
            environment( build_func ) <- parent.env(environment())
            out <- do.call( build_func , 
                    list( left , right , dist_name="uniform" , suffix="lpdf" , as_log_lik=as_log_lik , vectorized=TRUE ) , 
                    envir=environment() , quote=TRUE )

            return( out )
        }
    ) ,
    normal = list(
        R_name = "dnorm",
        Stan_name = "normal",
        Stan_suffix = "lpdf",
        pars = 2,
        dims = c( "real" , "real" , "real" ),
        constraints = c( NA , NA , "lower=0" ),
        vectorized = TRUE,
        build = function( left , right , as_log_lik=FALSE ) {
            #right <- as.character( right )
            #out <- ""

            # check for matrix var as outcome
            # in that case, use to_vector to vectorize assignment
            if ( length(left)==1 ) { # only applies to single symbol on left
                if ( !is.null(symbols[[left]]$dims) ) {
                    dims <- symbols[[left]]$dims
                    if ( length(dims) > 1 ) {
                        # list, so check if matrix with dimensions
                        if ( dims[[1]]=="matrix" ) {
                            left <- concat( "to_vector( " , left , " )" )
                        }
                        # check if array --- currently arrays have character dims like '2,2,2'
                        if ( dims[[1]]=="real" | dims[[1]]=="int" ) {
                            if ( class(dims[[2]])=="character" ) {
                                # see if contains a comma
                                sr <- strsplit( dims[[2]] , "," )[[1]]
                                if ( length(sr)>1 ) {
                                    # contains at least one comma
                                    # so need to loop over dimensions to assign prior
                                    n_idx <- length(sr)
                                    sr <- as.numeric(sr) # hope this works!
                                    loop_txt <- concat( "for ( ix1 in 1:" , sr[1] , " ) " )
                                    ix_txt <- "ix1"
                                    if ( n_idx > 1 ) {
                                        for ( xx in 2:n_idx ) {
                                            loop_txt <- concat( loop_txt , "for ( ix",xx," in 1:" , sr[xx] , " ) " )
                                            ix_txt <- concat( ix_txt , "," , "ix" , xx )
                                        }#xx
                                    }
                                    left <- concat( loop_txt , left , "[" , ix_txt , "]" )
                                }
                            }
                        }
                    }#length(dims)>1
                }
            }

            #for ( j in 1:length(left) )
            #    out <- concat( out , "    " , left[j] , " ~ " , "normal" , "( " , right[2] , " , " , right[3] , " );\n" )

            # call DEFAULT function for rest
            #out <- distribution_library[["DEFAULT"]]$build( left , right , dist_name="normal" , suffix="lpdf" , as_log_lik=as_log_lik )

            build_func <- distribution_library[["DEFAULT"]]$build
            environment( build_func ) <- parent.env(environment())
            out <- do.call( build_func , 
                    list( left , right , dist_name="normal" , suffix="lpdf" , as_log_lik=as_log_lik , vectorized=TRUE ) , 
                    envir=environment() , quote=TRUE )

            return( out )
        }
    ) ,
    # half_normal is just normal with a <lower=0> constraint on outcome
    half_normal = list(
        R_name = "dhalfnorm",
        Stan_name = "normal",
        Stan_suffix = "lpdf",
        pars = 2,
        dims = c( "real" , "real" , "real" ),
        constraints = c( "lower=0" , NA , "lower=0" ),
        vectorized = TRUE
    ) ,
    log_normal = list(
        R_name = "dlnorm",
        Stan_name = "lognormal",
        Stan_suffix = "lpdf",
        pars = 2,
        dims = c( "real" , "real" , "real" ),
        constraints = c( "lower=0" , NA , "lower=0" ),
        vectorized = TRUE
    ) ,
    # multi_normal can take (mu,Rho,sigma) or (mu,Sigma) or (mu,L_Sigma) as arguments
    # i.e. is overloaded
    multinormal1 = list(
        R_name = "dmvnorm",
        Stan_name = "multi_normal",
        Stan_suffix = "lpdf",
        pars = 2,
        dims = c( "vector" , "vector" , "matrix" , "vector" ),
        constraints = c( NA , NA , NA , "lower=0" ),
        vectorized = TRUE,
        build = function( left , right , as_log_lik=FALSE ) {

            flag_data_outcome <- FALSE
            if ( left[1] %in% names(data) ) flag_data_outcome <- TRUE

            # need to figure out dimension of distribution
            # can do so from either length of vector on left or length of MU or size of SIGMA/RHO
            n_vars <- length(left)
            n_cases <- 1
            if ( n_vars > 1 ) {
                # two possibilities:
                # (1) dims could have any grouping variable to tell us how many cases
                # (2) dims are just the rows in the table
                if ( !is.null(attr(left,"dims")) ) {
                    # must be a grouping variable?
                    n_cases <- length( unique( data[[ as.character( attr(left,"dims") ) ]] ) )
                } else {
                    # consider dims in symbol list
                    the_dims <- symbols[[ left[1] ]]$dims
                    # fetch case count
                    n_cases <- as.numeric( the_dims[[2]] )
                }
            }
            
            SIGMA <- as.character( right[[3]] )
            if ( length(right)==4 ) {
                # correlation matrix and vector of scales
                # so compose into covariance matrix
                SIGMA <- concat( "quad_form_diag(" , as.character(right[[3]]) , " , " , as.character(right[[4]]) , ")" )
            }

            MU <- as.character( right[[2]] )
            # check for c() and cbind() style vector of means
            if ( class( right[[2]] )=="call" ) {
                if ( as.character(right[[2]][[1]]) %in% c( "c" , "cbind" ) ) {
                    MU <- MU[-1] # -1 removes "c"
                    n_vars <- length(MU)
                } else {
                    stop( concat( "Something unexpected: " , deparse(right) ) )
                }
            }

            if ( n_vars==1 ) {
                # still haven't figured out dimensions
                # check explicit dims on left side variable
                if ( !is.null(attr(left,"dims")) ) {
                    dims <- attr(left,"dims")
                    if ( dims[[1]]=="vector" )
                        n_vars <- dims[[2]]
                    # check cases too
                    if ( length(dims)==3 ) {
                        n_cases <- length( unique( data[[ as.character( dims[[3]] ) ]] ) )
                    }
                }
            }

            #print(n_vars)
            #print(n_cases)

            # try to set dimensions on any Rho, sigma, SIGMA parameters
            if ( length(right)==4 ) {
                # have Rho, sigma as parameters
                Rho_name <- as.character( right[[3]] )
                sigma_name <- as.character( right[[4]] )
                # try to find in symbols list (in parent calling environment)
                # should have access, since build() is called with parent.env(environment())
                if ( !is.null( symbols[[Rho_name]] ) ) {
                    dims <- symbols[[Rho_name]]$dims
                    if ( !is.null(dims) ) {
                        if ( class(dims)=="character" || length(dims)==1 ) {
                            new_symbols <- symbols
                            new_symbols[[Rho_name]]$dims <- list( "corr_matrix" , n_vars )
                            assign( "symbols" , new_symbols , inherits=TRUE )
                        }
                    } else {
                        # null so do all the dims
                        new_symbols <- symbols
                        new_symbols[[Rho_name]]$dims <- list( "corr_matrix" , n_vars )
                        assign( "symbols" , new_symbols , inherits=TRUE )
                    }
                }
                if ( !is.null( symbols[[sigma_name]] ) ) {
                    dims <- symbols[[sigma_name]]$dims
                    if ( !is.null(dims) ) {
                        if ( class(dims)=="character" || length(dims)==1 ) {
                            new_symbols <- symbols
                            new_symbols[[sigma_name]]$dims <- list( "vector" , n_vars )
                            assign( "symbols" , new_symbols , inherits=TRUE )
                        }
                    } else {
                        # null so do all the dims
                        new_symbols <- symbols
                        new_symbols[[sigma_name]]$dims <- list( "vector" , n_vars )
                        assign( "symbols" , new_symbols , inherits=TRUE )
                    }
                }
            } else {
                # have only SIGMA covariance parameter
                SIGMA_name <- as.character( right[[3]] )
                if ( !is.null( symbols[[SIGMA_name]] ) ) {
                    dims <- symbols[[SIGMA_name]]$dims
                    if ( !is.null(dims) ) {
                        if ( class(dims)=="character" || length(dims)==1 ) {
                            new_symbols <- symbols
                            new_symbols[[SIGMA_name]]$dims <- list( "cov_matrix" , n_vars )
                            assign( "symbols" , new_symbols , inherits=TRUE )
                        }
                    } else {
                        # null so do all the dims
                        new_symbols <- symbols
                        new_symbols[[SIGMA_name]]$dims <- list( "cov_matrix" , n_vars )
                        assign( "symbols" , new_symbols , inherits=TRUE )
                    }
                }
            }

            # check for zero mean
            if ( length(MU)==1 )
                if ( MU=="0" ) {
                    MU <- concat( "rep_vector(0," , n_vars , ")" )
                }

            # now build text
            # local model block code constructs vector for left side
            #{
            #    vector[2] YY[6];
            #    vector[2] MU;
            #    MU = [mu1,mu2]';
            #    for ( j in 1:6 )
            #        YY[j] = [ a[j] , b[j] ]';
            #    YY ~ multi_normal_lpdf( MU , quad_form_diag(Rho,sigma) );
            #}

            out <- ""
            indent <- "    "
            # do we need a {} environment?
            if ( length(MU)>1 || length(left)>1 ) {
                out <- concat( out , indent , "{\n" )
            }

            # do we need a local var for left side?
            if ( length(left) > 1 ) {
                out <- concat( out , indent , "vector[" , n_vars , "] YY" )
                if ( n_cases > 1 )
                    out <- concat( out , "[" , n_cases , "];\n" )
                else
                    out <- concat( out , ";\n" )
            }

            # do we need a local var for means?
            if ( length(MU)>1 ) {
                # what is length of each element?
                # should either be 1 or same as n_cases
                vlen <- 1
                vdims <- symbols[[ MU[1] ]]$dims
                if ( class(vdims)=="list" ) {
                    vlen <- vdims[[2]]
                    if ( vlen != n_cases ) warning( "multi_normal mean vector has length > 1 but not same length as outcome" )
                }
                # build text
                out <- concat( out , indent , "vector[" , n_vars , "] MU" )
                if ( vlen > 1 ) 
                    out <- concat( out , "[" , vlen , "];\n" )
                else
                    out <- concat( out , ";\n" )
                # assign it too
                vsuf <- ""
                if ( vlen==1 )
                    out <- concat( out , indent , "MU = [ " )
                else {
                    out <- concat( out , indent , "for ( j in 1:" , vlen , " ) " )
                    out <- concat( out , "MU[j] = [ " )
                    vsuf <- "[j]"
                }
                for ( j in 1:length(MU) ) {
                    out <- concat( out , MU[j] , vsuf )
                    if ( j < length(MU) ) out <- concat( out , " , " )
                }
                out <- concat( out , " ]';\n" )
            }

            # loop to populate local for left side
            if ( length(left)>1 && n_cases>1 ) {
                out <- concat( out , indent , "for ( j in 1:" , n_cases , " ) " )
                out <- concat( out , "YY[j] = [ " )
                for ( j in 1:n_vars ) {
                    out <- concat( out , left[j] , "[j]" )
                    if ( j < n_vars ) out <- concat( out , " , " )
                }
                out <- concat( out , " ]';\n")
            }

            # finally, the distribution statement
            
            out_var <- "YY"
            if ( length(left)==1 ) out_var <- left[1]
            l1dims <- symbols[[ left[1] ]]$dims
            out_mat <- FALSE
            if ( class(l1dims) != "name" ) {
                if ( l1dims[[1]]=="matrix" ) out_mat <- TRUE
            }
            if ( flag_data_outcome==TRUE || out_mat==TRUE ) {
                if ( out_mat==TRUE ) {
                    N_out <- symbols[[ left[1] ]]$dims[[2]]
                    out_var <- concat( "for ( i in 1:" , N_out , " ) " , out_var , "[i,:]" )
                }
            }
            MU_var <- "MU"
            if ( length(MU)==1 ) MU_var <- MU[1]

            if ( as_log_lik==FALSE ) {

                out <- concat( out , indent , out_var , " ~ multi_normal( " , MU_var , " , " , SIGMA , " );\n" )

            } else {
                out <- concat( out , indent , "log_lik = multi_normal_lpdf( " , out_var , " | " , MU_var , " , " , SIGMA , " );\n" )
            }

            # close local environment, if necessary
            if ( length(MU)>1 || length(left)>1 ) {
                out <- concat( out , indent , "}\n" )
            }


            return( out )
        }
    ) ,
    # multi_normal that takes cholecky factor instead of covariance matrix
    multinormal_chol = list(
        R_name = "dmvnorm",
        Stan_name = "multi_normal_cholesky",
        Stan_suffix = "lpdf",
        pars = 2,
        dims = c( "vector" , "vector" , "matrix" ),
        constraints = c( NA , NA , NA ),
        vectorized = TRUE,
        build = function( left , right , as_log_lik=FALSE ) {

            flag_data_outcome <- FALSE
            if ( left[1] %in% names(data) ) flag_data_outcome <- TRUE

            # need to figure out dimension of distribution
            # can do so from either length of vector on left or length of MU or size of SIGMA/RHO
            n_vars <- length(left)
            n_cases <- 1
            if ( n_vars > 1 ) {
                # two possibilities:
                # (1) dims could have any grouping variable to tell us how many cases
                # (2) dims are just the rows in the table
                if ( !is.null(attr(left,"dims")) ) {
                    # must be a grouping variable?
                    n_cases <- length( unique( data[[ as.character( attr(left,"dims") ) ]] ) )
                } else {
                    # consider dims in symbol list
                    the_dims <- symbols[[ left[1] ]]$dims
                    # fetch case count
                    n_cases <- as.numeric( the_dims[[2]] )
                }
            }
            
            SIGMA <- as.character( right[[3]] )
            if ( length(right)==4 ) {
                # correlation matrix and vector of scales
                # so compose into covariance matrix
                SIGMA <- concat( "quad_form_diag(" , as.character(right[[3]]) , " , " , as.character(right[[4]]) , ")" )
            }

            MU <- as.character( right[[2]] )
            # check for c() and cbind() style vector of means
            if ( class( right[[2]] )=="call" ) {
                if ( as.character(right[[2]][[1]]) %in% c( "c" , "cbind" ) ) {
                    MU <- MU[-1] # -1 removes "c"
                    n_vars <- length(MU)
                } else {
                    stop( concat( "Something unexpected: " , deparse(right) ) )
                }
            }

            if ( n_vars==1 ) {
                # still haven't figured out dimensions
                # check explicit dims on left side variable
                if ( !is.null(attr(left,"dims")) ) {
                    dims <- attr(left,"dims")
                    if ( dims[[1]]=="vector" )
                        n_vars <- dims[[2]]
                    # check cases too
                    if ( length(dims)==3 ) {
                        n_cases <- length( unique( data[[ as.character( dims[[3]] ) ]] ) )
                    }
                }
            }

            #print(n_vars)
            #print(n_cases)

            # try to set dimensions on any Rho, sigma, SIGMA parameters
            if ( length(right)==4 ) {
                # have Rho, sigma as parameters
                Rho_name <- as.character( right[[3]] )
                sigma_name <- as.character( right[[4]] )
                # try to find in symbols list (in parent calling environment)
                # should have access, since build() is called with parent.env(environment())
                if ( !is.null( symbols[[Rho_name]] ) ) {
                    dims <- symbols[[Rho_name]]$dims
                    if ( !is.null(dims) ) {
                        if ( class(dims)=="character" || length(dims)==1 ) {
                            new_symbols <- symbols
                            new_symbols[[Rho_name]]$dims <- list( "corr_matrix" , n_vars )
                            assign( "symbols" , new_symbols , inherits=TRUE )
                        }
                    } else {
                        # null so do all the dims
                        new_symbols <- symbols
                        new_symbols[[Rho_name]]$dims <- list( "corr_matrix" , n_vars )
                        assign( "symbols" , new_symbols , inherits=TRUE )
                    }
                }
                if ( !is.null( symbols[[sigma_name]] ) ) {
                    dims <- symbols[[sigma_name]]$dims
                    if ( !is.null(dims) ) {
                        if ( class(dims)=="character" || length(dims)==1 ) {
                            new_symbols <- symbols
                            new_symbols[[sigma_name]]$dims <- list( "vector" , n_vars )
                            assign( "symbols" , new_symbols , inherits=TRUE )
                        }
                    } else {
                        # null so do all the dims
                        new_symbols <- symbols
                        new_symbols[[sigma_name]]$dims <- list( "vector" , n_vars )
                        assign( "symbols" , new_symbols , inherits=TRUE )
                    }
                }
            } else {
                # have only SIGMA covariance parameter
                SIGMA_name <- as.character( right[[3]] )
                if ( !is.null( symbols[[SIGMA_name]] ) ) {
                    dims <- symbols[[SIGMA_name]]$dims
                    if ( !is.null(dims) ) {
                        if ( class(dims)=="character" || length(dims)==1 ) {
                            new_symbols <- symbols
                            new_symbols[[SIGMA_name]]$dims <- list( "cov_matrix" , n_vars )
                            assign( "symbols" , new_symbols , inherits=TRUE )
                        }
                    } else {
                        # null so do all the dims
                        new_symbols <- symbols
                        new_symbols[[SIGMA_name]]$dims <- list( "cov_matrix" , n_vars )
                        assign( "symbols" , new_symbols , inherits=TRUE )
                    }
                }
            }

            # check for zero mean
            if ( length(MU)==1 )
                if ( MU=="0" ) {
                    MU <- concat( "rep_vector(0," , n_vars , ")" )
                }

            # now build text
            # local model block code constructs vector for left side
            #{
            #    vector[2] YY[6];
            #    vector[2] MU;
            #    MU = [mu1,mu2]';
            #    for ( j in 1:6 )
            #        YY[j] = [ a[j] , b[j] ]';
            #    YY ~ multi_normal_lpdf( MU , quad_form_diag(Rho,sigma) );
            #}

            out <- ""
            indent <- "    "
            # do we need a {} environment?
            if ( length(MU)>1 || length(left)>1 ) {
                out <- concat( out , indent , "{\n" )
            }

            # do we need a local var for left side?
            if ( length(left) > 1 ) {
                out <- concat( out , indent , "vector[" , n_vars , "] YY" )
                if ( n_cases > 1 )
                    out <- concat( out , "[" , n_cases , "];\n" )
                else
                    out <- concat( out , ";\n" )
            }

            # do we need a local var for means?
            if ( length(MU)>1 ) {
                # what is length of each element?
                # should either be 1 or same as n_cases
                vlen <- 1
                vdims <- symbols[[ MU[1] ]]$dims
                if ( class(vdims)=="list" ) {
                    vlen <- vdims[[2]]
                    if ( vlen != n_cases ) warning( "multi_normal mean vector has length > 1 but not same length as outcome" )
                }
                # build text
                out <- concat( out , indent , "vector[" , n_vars , "] MU" )
                if ( vlen > 1 ) 
                    out <- concat( out , "[" , vlen , "];\n" )
                else
                    out <- concat( out , ";\n" )
                # assign it too
                vsuf <- ""
                if ( vlen==1 )
                    out <- concat( out , indent , "MU = [ " )
                else {
                    out <- concat( out , indent , "for ( j in 1:" , vlen , " ) " )
                    out <- concat( out , "MU[j] = [ " )
                    vsuf <- "[j]"
                }
                for ( j in 1:length(MU) ) {
                    out <- concat( out , MU[j] , vsuf )
                    if ( j < length(MU) ) out <- concat( out , " , " )
                }
                out <- concat( out , " ]';\n" )
            }

            # loop to populate local for left side
            if ( length(left)>1 && n_cases>1 ) {
                out <- concat( out , indent , "for ( j in 1:" , n_cases , " ) " )
                out <- concat( out , "YY[j] = [ " )
                for ( j in 1:n_vars ) {
                    out <- concat( out , left[j] , "[j]" )
                    if ( j < n_vars ) out <- concat( out , " , " )
                }
                out <- concat( out , " ]';\n")
            }

            # finally, the distribution statement
            
            out_var <- "YY"
            if ( length(left)==1 ) out_var <- left[1]
            l1dims <- symbols[[ left[1] ]]$dims
            out_mat <- FALSE
            if ( class(l1dims) != "name" ) {
                if ( l1dims[[1]]=="matrix" ) out_mat <- TRUE
            }
            if ( flag_data_outcome==TRUE || out_mat==TRUE ) {
                if ( out_mat==TRUE ) {
                    N_out <- symbols[[ left[1] ]]$dims[[2]]
                    out_var <- concat( "for ( i in 1:" , N_out , " ) " , out_var , "[i,:]" )
                }
            }
            MU_var <- "MU"
            if ( length(MU)==1 ) MU_var <- MU[1]

            if ( as_log_lik==FALSE ) {

                out <- concat( out , indent , out_var , " ~ multi_normal_cholesky( " , MU_var , " , " , SIGMA , " );\n" )

            } else {
                out <- concat( out , indent , "log_lik = multi_normal_cholesky_lpdf( " , out_var , " | " , MU_var , " , " , SIGMA , " );\n" )
            }

            # close local environment, if necessary
            if ( length(MU)>1 || length(left)>1 ) {
                out <- concat( out , indent , "}\n" )
            }


            return( out )
        }
    ) ,
    bernoulli = list(
        R_name = "dbern",
        Stan_name = "bernoulli",
        Stan_suffix = "lpmf",
        pars = 1,
        dims = c( "int" , "real" ),
        constraints = c( NA , "lower=0,upper=1" ),
        vectorized = TRUE
    ),
    bernoulli_logit = list(
        R_name = "dbern",
        Stan_name = "bernoulli_logit",
        Stan_suffix = "lpmf",
        pars = 1,
        dims = c( "int" , "real" ),
        constraints = c( NA , "lower=0,upper=1" ),
        vectorized = TRUE
    ),
    binomial = list(
        R_name = "dbinom",
        Stan_name = "binomial",
        Stan_suffix = "lpmf",
        pars = 2,
        dims = c( "int" , "int" , "real" ),
        constraints = c( NA , NA , "lower=0,upper=1" ),
        vectorized = TRUE
    ),
    betabinomial = list(
        # Stan uses alpha,beta parameterization
        # need to make this into p,theta parameterization
        # p = alpha/(alpha+beta)
        # theta = alpha + beta
        R_name = "dbetabinom",
        Stan_name = "beta_binomial",
        Stan_suffix = "lpmf",
        pars = 3,
        dims = c( "int" , "int" , "real" , "real" ),
        constraints = c( NA , "lower=0" , "lower=0" ),
        vectorized = TRUE ,
        build = function( left , right , as_log_lik=FALSE ) {

            # convert between parameterizations
            # alpha = p*theta, beta = (1-p)*theta

            # right: [[1]] dist, [[2]] size, [[3]] p, [[4]] theta
            p_name <- as.character( right[[3]] )
            theta_name <- as.character( right[[4]] )
            right[[3]] <- concat( p_name , "*" , theta_name )
            right[[4]] <- concat( "(1-" , p_name , ")*" , theta_name )
            
            # call DEFAULT function for rest
            build_func <- distribution_library[["DEFAULT"]]$build
            environment( build_func ) <- parent.env(environment())
            out <- do.call( build_func , 
                    list( left , right , dist_name="beta_binomial" , suffix="lpmf" , as_log_lik=as_log_lik , vectorized=TRUE ) , 
                    envir=environment() , quote=TRUE )

            return( out )
        }
    ),
    gamma2 = list(
        # Stan uses shape,inv-scale parameterization
        # need to make this into mu,scale parameterization
        # alpha = mu/theta
        # beta = 1/theta
        R_name = "dgamma2",
        Stan_name = "gamma",
        Stan_suffix = "lpdf",
        pars = 2,
        dims = c( "real" , "real" , "real" ),
        constraints = c( "lower=0" , "lower=0" , "lower=0" ),
        vectorized = TRUE ,
        build = function( left , right , as_log_lik=FALSE ) {

            # convert between parameterizations

            # right: [[1]] dist, [[2]] mu [[3]] theta
            mu_name <- as.character( right[[2]] )
            theta_name <- as.character( right[[3]] )
            right[[2]] <- concat( mu_name , "/" , theta_name )
            right[[3]] <- concat( "1/" , theta_name )
            
            # call DEFAULT function for rest
            build_func <- distribution_library[["DEFAULT"]]$build
            environment( build_func ) <- parent.env(environment())
            out <- do.call( build_func , 
                    list( left , right , dist_name="gamma" , suffix="lpdf" , as_log_lik=as_log_lik , vectorized=TRUE ) , 
                    envir=environment() , quote=TRUE )

            return( out )
        }
    ),
    binomial_logit = list(
        R_name = "dbinom",
        Stan_name = "binomial_logit",
        Stan_suffix = "lpmf",
        pars = 2,
        dims = c( "int" , "int" , "real" ),
        constraints = c( NA , NA , NA ),
        vectorized = TRUE
    ),
    poisson = list(
        R_name = "dpois",
        Stan_name = "poisson",
        Stan_suffix = "lpmf",
        pars = 1,
        dims = c( "int" , "real" ),
        constraints = c( NA , "lower=0" ),
        vectorized = TRUE
    ),
    gammapoisson = list(
        R_name = "dgampois",
        Stan_name = "neg_binomial_2", # mu/phi parameterization
        Stan_suffix = "lpmf",
        pars = 2,
        dims = c( "int" , "real" , "real" ),
        constraints = c( NA , "lower=0" , "lower=0" ),
        vectorized = TRUE
    ),
    exponential = list(
        R_name = "dexp",
        Stan_name = "exponential",
        Stan_suffix = "lpdf",
        pars = 1,
        dims = c( "real" , "real" ),
        constraints = c( "lower=0" , "lower=0" ),
        vectorized = TRUE
    ),
    pareto1 = list(
        R_name = "dpareto",
        Stan_name = "pareto",
        Stan_suffix = "lpdf",
        pars = 2,
        dims = c( "real" , "real" ),
        constraints = c( "lower=0" , "lower=0" ),
        vectorized = TRUE
    ),
    beta = list(
        R_name = "dbeta",
        Stan_name = "beta",
        Stan_suffix = "lpdf",
        pars = 2,
        dims = c( "real" , "real" , "real" ),
        constraints = c( "lower=0,upper=1" , "lower=0" , "lower=0" ),
        vectorized = TRUE
    ),
    # version of beta with p,theta
    # p = alpha/(alpha+beta)
    # theta = alpha + beta
    beta2 = list(
        R_name = "dbeta2",
        Stan_name = "beta",
        Stan_suffix = "lpdf",
        pars = 2,
        dims = c( "real" , "real" , "real" ),
        constraints = c( "lower=0,upper=1" , "lower=0,upper=1" , "lower=0" ),
        vectorized = TRUE,
        build = function( left , right , as_log_lik=FALSE ) {

            # convert between parameterizations
            # alpha = p*theta, beta = (1-p)*theta

            # right: [[1]] dist, [[2]] p, [[3]] theta
            p_name <- as.character( right[[2]] )
            theta_name <- as.character( right[[3]] )
            right[[2]] <- concat( p_name , "*" , theta_name )
            right[[3]] <- concat( "(1-" , p_name , ")*" , theta_name )
            
            # call DEFAULT function for rest
            build_func <- distribution_library[["DEFAULT"]]$build
            environment( build_func ) <- parent.env(environment())
            out <- do.call( build_func , 
                    list( left , right , dist_name="beta" , suffix="lpdf" , as_log_lik=as_log_lik , vectorized=TRUE ) , 
                    envir=environment() , quote=TRUE )

            return( out )
        }
    ),
    dirichlet = list(
        R_name = "ddirichlet", # in gtools
        Stan_name = "dirichlet",
        Stan_suffix = "lpdf",
        pars = 1,
        dims = c( "simplex" , "vector" ),
        constraints = c( "lower=0,upper=1" , "lower=0" ),
        vectorized = TRUE
    ),
    cauchy = list(
        R_name = "dcauchy",
        Stan_name = "cauchy",
        Stan_suffix = "lpdf",
        pars = 2,
        dims = c( "real" , "real" , "real" ),
        constraints = c( NA , NA , "lower=0" ),
        vectorized = TRUE
    ),
    # half_cauchy is just cauchy with a <lower=0> constraint on outcome
    half_cauchy = list(
        R_name = "dhalfcauchy",
        Stan_name = "cauchy",
        Stan_suffix = "lpdf",
        pars = 2,
        dims = c( "real" , "real" , "real" ),
        constraints = c( "lower=0" , NA , "lower=0" ),
        vectorized = TRUE
    ) ,
    student_t = list(
        R_name = "dstudent",
        Stan_name = "student_t",
        Stan_suffix = "lpdf",
        pars = 3, # nu, mu, sigma
        dims = c( "real" , "real" , "real" , "real" ),
        constraints = c( NA , NA , NA , "lower=0" ),
        vectorized = TRUE
    ),
    ordered_logistic = list(
        R_name = "dordlogit",
        Stan_name = "ordered_logistic",
        Stan_suffix = "lpmf",
        pars = 2,
        dims = c( "int" , "real" , "ordered" ),
        constraints = c( "lower=1" , NA , NA ),
        vectorized = FALSE,
        build = function( left , right , as_log_lik=FALSE ) {

            # check type for vector of cutpoints
            # coerce to ordered of proper length, if no vector type already specified
            # cutpoints vector should be right[[3]]
            if ( !is.null( symbols[[ as.character(right[[3]]) ]] ) ) {
                cuts_name <- as.character(right[[3]])
                if ( class(symbols[[ cuts_name ]]$dims)=="character" ) {
                    # prob defaulted to "real" - set to ordered type
                    n_levels <- max( data[[ left ]] ) # can only be data
                    new_symbols <- symbols
                    new_symbols[[ cuts_name ]]$dims <- list( "ordered" , n_levels-1 )
                    assign( "symbols" , new_symbols , inherits=TRUE )
                }
            }
            
            # call DEFAULT function for rest
            #out <- distribution_library[["DEFAULT"]]$build( left , right , dist_name="ordered_logistic" , suffix="lpmf" , as_log_lik=as_log_lik , vectorized=FALSE )
            build_func <- distribution_library[["DEFAULT"]]$build
            environment( build_func ) <- parent.env(environment())
            out <- do.call( build_func , 
                    list( left , right , dist_name="ordered_logistic" , suffix="lpmf" , as_log_lik=as_log_lik , vectorized=FALSE ) , 
                    envir=environment() , quote=TRUE )

            return( out )
        }
    ),
    categorical = list(
        R_name = "dcategorical",
        Stan_name = "categorical",
        Stan_suffix = "lpmf",
        pars = 1,
        dims = c( "int" , "vector" ),
        constraints = c( NA , NA ),
        vectorized = TRUE
    ),
    categorical_logit = list(
        R_name = "dcategorical",
        Stan_name = "categorical_logit",
        Stan_suffix = "lpmf",
        pars = 1,
        dims = c( "int" , "vector" ),
        constraints = c( NA , NA ),
        vectorized = TRUE
    ),
    lkjcorr = list(
        R_name = "dlkjcorr",
        Stan_name = "lkj_corr",
        Stan_suffix = "lpdf",
        pars = 1,
        dims = c( "corr_matrix" , "real" ),
        constraints = c( NA , "lower=0" ),
        vectorized = TRUE
    ),
    lkjcorrchol = list(
        R_name = "dlkjcorr",
        Stan_name = "lkj_corr_cholesky",
        Stan_suffix = "lpdf",
        pars = 1,
        dims = c( "cholesky_factor_corr" , "real" ),
        constraints = c( NA , "lower=0" ),
        vectorized = TRUE
    ),
    zipoisson = list(
        R_name = "dzipois",
        Stan_name = "poisson",
        Stan_suffix = "lpmf",
        pars = 2,
        dims = c( "int" , "real" , "real" ),
        constraints = c( "lower=0" , "lower=0,upper=1" , "lower=0" ),
        vectorized = FALSE,
        build = function( left , right , as_log_lik=FALSE ) {
            # need to build code for multiple choice likelihood
            # Pr(y|y==0) = p + (1-p)*Pr(y==0|lambda)
            # Pr(y|y>0)  = (1-p)*Pr(y|lambda)
            # where p is probability of zero inflation process (can be a local variable)
            # and lambda is poisson rate (can be local variable)

            right <- as.character( right )
            out <- ""
            p_symbol <- right[2]
            lambda_symbol <- right[3]
            if ( !is.null(symbols[[ p_symbol ]]) ) {
                if ( symbols[[ p_symbol ]]$type %in% c("local","data") ) {
                    if ( class(symbols[[ p_symbol ]]$dims)=="list" )
                        p_symbol <- concat( p_symbol , "[i]" )
                }
            }
            if ( !is.null(symbols[[ lambda_symbol ]]) ) {
                if ( symbols[[ lambda_symbol ]]$type %in% c("local","data") ) {
                    if ( class(symbols[[ lambda_symbol ]]$dims)=="list" )
                        lambda_symbol <- concat( lambda_symbol , "[i]" )
                }
            }

            out_head <- "target += "
            if ( as_log_lik==TRUE ) out_head <- "log_lik[i] = "

            out <- concat( out , indent , "for ( i in 1:" , length(data[[ left ]]) , " ) {\n" )

            out <- concat( out , indent , indent , "if ( " , left , "[i]==0 )\n" )
            out <- concat( out , inden(3) , out_head , "log_mix( " , p_symbol , " , 0 , poisson_lpmf(0|" , lambda_symbol , ") );\n" )
            out <- concat( out , indent , indent , "if ( " , left , "[i] > 0 )\n" )
            out <- concat( out , inden(3) , out_head , "log1m( " , p_symbol , " ) + poisson_lpmf(" , left , "[i] | " , lambda_symbol , " );\n" )

            out <- concat( out , indent , "}\n" )

            return(out)
        }
    ),
    CJS_capture_history = list(
        R_name = "dCJScapture",
        Stan_name = "CJS_capture_history",
        Stan_suffix = "lpmf",
        pars = 2,
        dims = c( "int" , "vector" , "vector" ),
        constraints = c( "lower=0,upper=1" , "lower=0,upper=1" , "lower=0,upper=1" ),
        vectorized = FALSE,
        build = function( left , right , as_log_lik=FALSE ) {
            # need to gsub left-hand variable into the transformed data code 
            new_tdata1 <- gsub( "CJS_OUTCOME" , left , m_tdata1 , fixed=TRUE )
            assign( "m_tdata1" , new_tdata1 , inherits=TRUE )
            new_tdata2 <- gsub( "CJS_OUTCOME" , left , m_tdata2 , fixed=TRUE )
            assign( "m_tdata2" , new_tdata2 , inherits=TRUE )

            # make sure outcome variable is typed as int[I,T], *not* as matrix
            # can be a matrix on R side --- Stan will coerce it to integer array
            old_dims <- symbols[[ left ]]$dims
            new_dims <- list( "int" , concat( old_dims[[2]] , "," , old_dims[[3]] ) )
            new_symbols <- symbols
            new_symbols[[ left ]]$dims <- new_dims
            assign( "symbols" , new_symbols , inherits=TRUE )

            # call DEFAULT function for rest
            build_func <- distribution_library[["DEFAULT"]]$build
            environment( build_func ) <- parent.env(environment())
            out <- do.call( build_func , 
                    list( left , right , dist_name="CJS_capture_history" , suffix="lpmf" , as_log_lik=as_log_lik , vectorized=FALSE ) , 
                    envir=environment() , quote=TRUE )
            # done
            return( out )
        }
    ),
    custom = list(
        # use '!Y' in formula to indicate location of outcome
        # or can hard code location of outcome
        R_name = "dcustom",
        Stan_name = "custom",
        Stan_suffix = "lpdf",
        pars = 1,
        dims = c( NA , NA ),
        constraints = c( NA , NA ),
        vectorized = FALSE,
        build = function( left , right , as_log_lik=FALSE ) {

            out <- ""
            vectorized <- FALSE

            for ( j in 1:length(left) ) {

                # check for any multiple choice flags on left symbol
                do_iff <- FALSE
                if( !is.null( attr( left , "iff" ) ) ) {
                    # need to add index to symbols in condition
                    iff_out <- attr( left , "iff" )
                    symlist <- unique( c( left[j] , names(data) , names(symbols) ) )
                    iff_out <- nest_indexes( iff_out , symlist )
                    do_iff <- TRUE
                }

                # in case of explicit loop, need max index
                # for matrices & arrays, loop over first index only
                max_index <- dim( data[[ left[j] ]] )[1]
                if ( is.null(max_index) ) max_index <- length( data[[ left[j] ]] )

                # need right side without the custom() call
                right_out <- right[[2]] # just inside custom()
                # check for '!Y'
                right_out <- hunt_bangY( right_out , left[j] )
                # first check symbols and add [i] where dim matches loop length
                sym <- get_all_symbols( right_out )
                if ( length(sym)>0 )
                    for( k in sym ) {
                        if ( k == left[j] ) 
                            right_out <- nest_indexes( right_out , left[j] )
                        else {
                            # not left symbol, so check dims
                            if ( k %in% names(symbols) ) {
                                if ( class(symbols[[k]]$dims) == "list" ) {
                                    # need length
                                    the_dim <- symbols[[k]]$dims[[2]]
                                    if ( the_dim == max_index )
                                        right_out <- nest_indexes( right_out , k )
                                }
                            }
                        }
                    }#k
                right_out <- deparse( right_out , width.cutoff=500 )
                if ( length(right_out)>1 ) right_out <- paste( right_out , collapse="" )

                # outcome var and distribution name
                if ( as_log_lik==FALSE ) {
                    out <- concat( out , "    " , "for ( i in 1:" , max_index , " ) " )
                    if ( do_iff == TRUE )
                        out <- concat( out , "\n" , inden(2) , "if ( " , deparse(iff_out) , " ) " )
                    out <- concat( out , "target += " , right_out , ";\n" )
                } else {
                    # log_lik expressions never vectorized, because vectorized distributions return a sum of log_likelihoods, and we need individual terms
                    out <- concat( out , "    " , "for ( i in 1:" , max_index , " ) " )
                    if ( do_iff == TRUE )
                        out <- concat( out , "\n" , inden(2) , "if ( " , deparse(iff_out) , " ) " )
                    out <- concat( out , "log_lik[i] = " , right_out , ";\n" )
                }

            }#j
            return( out )
        }
    ),
    mixture = list(
        # custom finite mixtures
        # format: y ~ mixture( theta , term1 , term2 , ... )
        # where theta has length 1 (for 2 terms) or same length as term count
        R_name = "dmixture",
        Stan_name = "mixture",
        Stan_suffix = "lpdf",
        pars = 2,
        dims = c( NA , NA , NA ),
        constraints = c( NA , NA , NA ),
        vectorized = FALSE,
        build = function( left , right , as_log_lik=FALSE ) {

            out <- ""
            vectorized <- FALSE

            for ( j in 1:length(left) ) {

                # check for any multiple choice flags on left symbol
                do_iff <- FALSE
                if( !is.null( attr( left , "iff" ) ) ) {
                    # need to add index to symbols in condition
                    iff_out <- attr( left , "iff" )
                    symlist <- unique( c( left[j] , names(data) , names(symbols) ) )
                    iff_out <- nest_indexes( iff_out , symlist )
                    do_iff <- TRUE
                }

                # in case of explicit loop, need max index
                # for matrices & arrays, loop over first index only
                max_index <- dim( data[[ left[j] ]] )[1]
                if ( is.null(max_index) ) max_index <- length( data[[ left[j] ]] )

                # need right side without the mixture() call

                # figure out how many terms we are mixing over
                theta <- right[[2]]
                n_terms <- length(right) - 2 # minus mixture and theta argument
                terms <- list()
                for ( i in 1:n_terms ) terms[[i]] <- right[[2+i]]

                # check for '!Y'
                #right_out <- hunt_bangY( right_out , left[j] )

                # process each term
                # all of these end up inside log_mix or log_sum_exp
                # so need to be ordinary log-prob terms to add to target
                for ( i in 1:n_terms ) {
                    # need to add _lpdf or _lpmf to any distribution names
                    front <- terms[[i]][[1]]
                    front_dist <- get_dist_template( as.character(front) )
                    terms[[i]][[1]] <- as.name( concat( front_dist$Stan_name , "_" , front_dist$Stan_suffix ) )

                    # also need to insert left symbol as first argument
                    # also need |
                    # can replace 2nd item term with `|` call
                    terms[[i]][[2]] <- call( "|" , as.name( left[j] ) , terms[[i]][[2]] )

                    # check symbols and add [i] where dim matches loop length
                    sym <- get_all_symbols( terms[[i]] )
                    if ( length(sym)>0 )
                        for( k in sym ) {
                            if ( k == left[j] ) 
                                terms[[i]] <- nest_indexes( terms[[i]] , left[j] )
                            else {
                                # not left symbol, so check dims
                                if ( k %in% names(symbols) ) {
                                    if ( class(symbols[[k]]$dims) == "list" ) {
                                        # need length
                                        the_dim <- symbols[[k]]$dims[[2]]
                                        if ( the_dim == max_index )
                                            terms[[i]] <- nest_indexes( terms[[i]] , k )
                                    }
                                }# in symbols
                                if ( k %in% names(data) ) {
                                    # might be data, but not used already in formula
                                    the_dim <- length( data[[k]] )
                                    if ( the_dim == max_index )
                                        terms[[i]] <- nest_indexes( terms[[i]] , k )
                                }
                            }
                        }#k
                    terms[[i]] <- deparse( terms[[i]] , width.cutoff=500 )
                    if ( length(terms[[i]])>1 ) terms[[i]] <- paste( terms[[i]] , collapse="" )
                }

                # now need to put the terms together inside either log_mix (2 terms) or log_sum_exp
                # if log_sum_exp, need to add log(theta[n]) to each
                if ( n_terms>2 ) right_out <- "log_sum_exp( "
                else {
                    # log_mix
                    # theta might be observed (measurement error)
                    theta_suffix <- ""
                    if ( !is.null( data[[ deparse(theta) ]] ) ) {
                        if ( length(data[[ deparse(theta) ]]) == max_index )
                            theta_suffix <- "[i]"
                    }
                    right_out <- concat( "log_mix( " , deparse(theta) , theta_suffix , " , " )
                }
                for ( i in 1:n_terms ) {
                    if ( n_terms>2 ) {
                        # add log(theta[n]) to front of term
                        log_pre <- concat( "log( " , deparse(theta) , "[" , i , "] ) + " )
                        terms[[i]] <- concat( log_pre , terms[[i]] )
                    }
                    right_out <- concat( right_out , terms[[i]] )
                    if ( i < n_terms ) right_out <- concat( right_out , " , " )
                }
                right_out <- concat( right_out , " )" )

                # outcome var and distribution name
                if ( as_log_lik==FALSE ) {
                    out <- concat( out , "    " , "for ( i in 1:" , max_index , " ) " )
                    if ( do_iff == TRUE )
                        out <- concat( out , "\n" , inden(2) , "if ( " , deparse(iff_out) , " ) " )
                    out <- concat( out , "target += " , right_out , ";\n" )
                } else {
                    # log_lik expressions never vectorized, because vectorized distributions return a sum of log_likelihoods, and we need individual terms
                    out <- concat( out , "    " , "for ( i in 1:" , max_index , " ) " )
                    if ( do_iff == TRUE )
                        out <- concat( out , "\n" , inden(2) , "if ( " , deparse(iff_out) , " ) " )
                    out <- concat( out , "log_lik[i] = " , right_out , ";\n" )
                }

            }#j
            return( out )
        }
    )
)#end of distribution library

ulam_macros <- list(
    # shortcut for awkward t(diag_pre_multiply( sigma , L_Rho ) * z)
    compose_noncentered = list(
        build = function( f , n ) {
            # f is entire formula
            # n is line with macro in it
            # t(diag_pre_multiply( sigma , L_Rho ) * z)
            sigma_symbol <- deparse( f[[n]][[3]][[2]] )
            corr_symbol <- deparse( f[[n]][[3]][[3]] )
            z_symbol <- deparse( f[[n]][[3]][[4]] )
            new_right <- concat( "t(diag_pre_multiply( " , sigma_symbol , " , " , corr_symbol , ") * ", z_symbol , " )" )
            f[[n]][[3]] <- parse( text=new_right )[[1]]
            return( f )
        }
    ),
    # shortcut for multiply_lower_tri_self_transpose(L_Rho)
    Chol_to_Corr = list(
        build = function( f , n ) {
            # f is entire formula
            # n is line with macro in it
            L_symbol <- deparse( f[[n]][[3]][[2]] )
            new_right <- concat( "multiply_lower_tri_self_transpose( " , L_symbol , " )" )
            f[[n]][[3]] <- parse( text=new_right )[[1]]
            return( f )
        }
    ),
    multi_normal_NC = list(
        build = function( f , n ) {
            # f is entire formula
            # n is line with macro in it
            # t(diag_pre_multiply( sigma , L_Rho ) * z)
            sigma_symbol <- deparse( f[[n]][[3]][[2]] )
            corr_symbol <- deparse( f[[n]][[3]][[1]] )
            # replace multi_normal_NC line with matrix assignment
            new_right <- concat( "t(diag_pre_multiply( " , sigma_symbol , " , " , corr_symbol , ") * z)" )
            f[[n]][[3]] <- parse( text=new_right )[[1]]
            # replace Rho with L_Rho declare
            return( f )
        }
    ),
    # define vector of observed and missing values
    # function does work during run time
    # but we also have to set right dims on merge vector and recode NA to a special value
    merge_missing = list(
        build = function( f , n ) {
            # f is entire formula
            # n is line with macro in it
            ff <- f[[n]] # line with merge_missing call

            # (1) make sure left side symbol is declared as vector of same length as x_obs
            merge_symbol <- get_left_symbol( ff )
            xobs_symbol <- as.character( ff[[3]] )[2]
            #print(str(merge_symbol))
            if ( is.null(attr(merge_symbol,"dims")) ) {
                N <- length( data[[ xobs_symbol ]] )
                new_left <- concat( "vector[",N,"]: " , merge_symbol )
                f[[n]][[2]] <- parse( text=new_left )[[1]]
            } else {
                N <- attr(merge_symbol,"dims")[[2]] # should be length of vector
            }

            # (2) replace NAs with a unique flag value - flag not used, just good for bookkeeping
            # and then construct a vector of indexes with missing values
            # add that vector to data list
            NA_idx <- which( is.na( data[[ xobs_symbol ]] ) )
            flag <- (-9)
            while ( any( data[[ xobs_symbol ]][-NA_idx] == flag ) ) {
                flag <- flag*10 - 9 # add another 9 on it
            }
            new_dat <- data
            new_dat[[ xobs_symbol ]][ NA_idx ] <- flag
            miss_indexes_name <- concat( xobs_symbol , "_missidx" )
            # as.array conversion to fix case with a single missing value
            # otherwise doesn't declare int[] properly
            # integer(1d) array is later detected in ulam() and declared as type "int_array"
            new_dat[[ miss_indexes_name ]] <- as.array( NA_idx )
            assign( "data" , new_dat , inherits=TRUE )

            # (3) insert the miss_indexes symbol into the call as first argument (Stan will need it)
            xmiss_symbol <- as.character( ff[[3]] )[3]
            new_right <- concat( "merge_missing( " , miss_indexes_name , " , to_vector(" , xobs_symbol , ") , " , xmiss_symbol , " )" )
            f[[n]][[3]] <- parse( text=new_right )[[1]]

            # (4) make sure x_impute symbol gets declared with length = num NAs
            # can just insert into symbols list
            impute_list <- list( type="par" , dims=list( "vector" , length(NA_idx) ) )
            new_symbols <- symbols
            new_symbols[[ xmiss_symbol ]] <- impute_list
            assign( "symbols" , new_symbols , inherits=TRUE )

            # (5) finally, make <- into <<- so later parsing knows to vectorize this assignment
            f[[n]][[1]] <- as.name( "<<-" )
            # done
            return( f )
        },
        functions = "
    vector merge_missing( int[] miss_indexes , vector x_obs , vector x_miss ) {
        int N = dims(x_obs)[1];
        int N_miss = dims(x_miss)[1];
        vector[N] merged;
        merged = x_obs;
        for ( i in 1:N_miss )
            merged[ miss_indexes[i] ] = x_miss[i];
        return merged;
    }"
    ),
    cov_GPL2 = list(
        functions = "
    matrix cov_GPL2(matrix x, real sq_alpha, real sq_rho, real delta) {
        int N = dims(x)[1];
        matrix[N, N] K;
        for (i in 1:(N-1)) {
          K[i, i] = sq_alpha + delta;
          for (j in (i + 1):N) {
            K[i, j] = sq_alpha * exp(-sq_rho * square(x[i,j]) );
            K[j, i] = K[i, j];
          }
        }
        K[N, N] = sq_alpha + delta;
        return K;
    }"
    ),
    cov_GPL1 = list(
        functions = "
    matrix cov_GPL1(matrix x, real sq_alpha, real sq_rho, real delta) {
        int N = dims(x)[1];
        matrix[N, N] K;
        for (i in 1:(N-1)) {
          K[i, i] = sq_alpha + delta;
          for (j in (i + 1):N) {
            K[i, j] = sq_alpha * exp(-sq_rho * x[i,j] );
            K[j, i] = K[i, j];
          }
        }
        K[N, N] = sq_alpha + delta;
        return K;
    }"
    ),
    cov_Pagel = list(
        functions = "
    matrix cov_Pagel(matrix x, real lambda) {
        int N = dims(x)[1];
        matrix[N, N] K = x;
        for (i in 1:(N-1)) {
          for (j in (i + 1):N) {
            K[i, j] = lambda * x[i,j];
            K[j, i] = K[i, j];
          }
        }
        return K;
    }"
    ),
    CJS_capture_history = list(
        # following page 221 of Stan manual (Individual Cormack-Jolly-Seber Model) 
        functions = "
    int first_capture(int[] y_i) {
        for (k in 1:size(y_i))
            if (y_i[k]) return k;
        return 0;
    }

    int last_capture(int[] y_i) {
        for (k_rev in 0:(size(y_i) - 1)) {
            int k;
            k = size(y_i) - k_rev;
            if (y_i[k]) 
                return k; 
        }
        return 0; 
    }

    vector prob_uncaptured( int T, vector p, vector phi) {
        vector[T] chi;
        chi[T] = 1.0;
        for ( t in 1:(T-1) ) {
            int t_curr;
            int t_next;
            t_curr = T - t;
            t_next = t_curr + 1;
            chi[t_curr] = (1-phi[t_curr]) + phi[t_curr] * (1-p[t_next]) * chi[t_next];
        }
        return chi;
    }

    // use like: for ( i in 1:I ) y[i] ~ CJS_capture_history( p , phi );
    real CJS_capture_history_lpmf( int[] y , vector p , vector phi ) {
        vector[ size(y) ] chi;
        int first;
        int last;
        int T;
        real lp;
        T = size(y);
        chi = prob_uncaptured( T , p , phi );
        lp = 0;
        first = first_capture(y);
        last = last_capture(y);
        if (first > 0) {
            for ( t in (first+1):last ) {
                lp = lp + bernoulli_lpmf( 1 | phi[t-1] );
                lp = lp + bernoulli_lpmf( y[t] | p[t] );
            }
            lp = lp + bernoulli_lpmf( 1 | chi[last] );
        }
        return lp;
    }",
        transformed_data_declare = "
    int T_CJS; // number of captures
    int I_CJS; // number of individuals
    vector<lower=0,upper=dims( CJS_OUTCOME )[1]>[ dims( CJS_OUTCOME )[2] ] n_captured;",
        transformed_data_body = "
    T_CJS = dims( CJS_OUTCOME )[2];
    I_CJS = dims( CJS_OUTCOME )[1];
    n_captured = rep_vector( 0 , T_CJS );
    for ( t in 1:T_CJS )
    for ( i in 1:I_CJS )
        if ( CJS_OUTCOME[i, t] )
            n_captured[t] = n_captured[t] + 1;"
    )
) #end of macro library

# R versions of macro functions, sometimes needed for link
# trick is these functions need to vectorize over parameters to be useful in link

#matrix cov_GPL2(matrix x, real sq_alpha, real sq_rho, real delta) {
#        int N = dims(x)[1];
#        matrix[N, N] K;
#        for (i in 1:(N-1)) {
#          K[i, i] = sq_alpha + delta;
#          for (j in (i + 1):N) {
#            K[i, j] = sq_alpha * exp(-sq_rho * square(x[i,j]) );
#            K[j, i] = K[i, j];
#          }
#        }
#        K[N, N] = sq_alpha + delta;
#        return K;
#    }
cov_GPL2_scal <- function( x , sq_alpha , sq_rho , delta ) {
    N <- ncol(x)
    K <- matrix( NA , N , N )
    for ( i in 1:(N-1) ) {
        K[i,i] <- sq_alpha + delta
        for ( j in (i+1):N ) {
            K[i,j] <- sq_alpha * exp( -sq_rho * x[i,j]^2 )
            K[j,i] <- K[i,j]
        }
    }
    K[N,N] <- sq_alpha + delta
    return(K)
}
cov_GPL2 <- Vectorize( cov_GPL2_scal , c("sq_alpha","sq_rho","delta") , SIMPLIFY=FALSE )
# z <- cov_GPL2( y3 , post$etasq[1:10] , post$rhosq[1:10] , 0.01 )
