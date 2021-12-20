# map2stan 3.0
# redesign for more direct to Stan
# less reliance on GLM structure
# allow explicit declarations
# see tests in ulam-tests.R for examples

# file argument:
# if NULL ignore
# if string X,
#   if file X.rds exists in working directory, load instead of fitting
#   otherwise save result as X.rds to working directory

# May 2020: 
# Added cmdstan argument to use cmdstanr interface when installed
# Added threads argument for reduce_sum conversion of model block
# threads > 1 only works with cmdstan=TRUE at moment

ulam_options <- new.env(parent=emptyenv())
ulam_options$use_cmdstan <- FALSE
if ( require(cmdstanr) ) ulam_options$use_cmdstan <- TRUE
set_ulam_cmdstan <- function(x=TRUE) assign( "use_cmdstan" , x , env=ulam_options )

ulam <- function( flist , data , pars , pars_omit , start , chains=1 , cores=1 , iter=1000 , warmup , control=list(adapt_delta=0.95) , distribution_library=ulam_dists , macro_library=ulam_macros , custom , constraints , declare_all_data=TRUE , log_lik=FALSE , sample=TRUE , messages=TRUE , pre_scan_data=TRUE , coerce_int=TRUE , sample_prior=FALSE , file=NULL , cmdstan=ulam_options$use_cmdstan , threads=1 , grain=1 , cpp_options=list() , cpp_fast=FALSE , rstanout=TRUE , ... ) {

    if ( !is.null(file) ) {
        rds_file_name <- concat( file , ".rds" )
        if ( file.exists(rds_file_name) ) {
            # load it
            result <- readRDS( file=rds_file_name )
            return(result) # exit
        }
    }

    if ( !missing(data) )
        data <- as.list(data)

    if ( missing(warmup) ) warmup <- floor(iter/2)

    # check for previous fit passed instead of formula
    prev_stanfit <- FALSE
    if ( class(flist)=="ulam" ) {
        prev_stanfit <- TRUE
        prev_stanfit_object <- flist@stanfit
        if ( missing(data) ) data <- flist@data
        flist <- flist@formula
    }

    # check for quap/map formula
    if ( class(flist)=="map" ) {
        if ( missing(data) ) data <- flist@data
        flist <- flist@formula
    }

    if ( pre_scan_data==TRUE ) {
        # pre-scan for character variables and remove
        x_to_remove <- c()
        for ( i in 1:length(data) ) {
            if ( class(data[[i]])[1] =="character" ) {
                x_to_remove <- c( x_to_remove , names(data)[i] )
            }
            if ( class(data[[i]])[1] =="factor" ) {
                if ( coerce_int==FALSE )
                    x_to_remove <- c( x_to_remove , names(data)[i] )
                else {
                    # coerce factor to index variable
                    data[[i]] <- as.integer( data[[i]] )
                }
            }
        }#i
        if ( length(x_to_remove)>0 ) {
            if ( messages==TRUE ) {
                message( "Removing one or more character or factor variables:" )
                message( x_to_remove )
            }
            for ( i in 1:length(data) )
                if ( names(data)[i] %in% x_to_remove ) data[[i]] <- NULL
        }
        # pre-scan for index variables (integer) that are numeric by accident
        for ( i in 1:length(data) ) {
            if ( class(data[[i]])[1]!="character" ) {
                if ( all( as.integer(data[[i]])==data[[i]] , na.rm=TRUE ) ) {
                    #data[[i]] <- as.integer(data[[i]])
                    if ( class(data[[i]])[1]!="integer" & !inherits(data[[i]],"matrix") ) {
                        if ( coerce_int==TRUE ) {
                            data[[i]] <- as.integer(data[[i]])
                        }
                        if ( messages==TRUE & coerce_int==FALSE ) {
                            vn <- names(data)[i]
                            message( "Cautionary note:" )
                            message( concat( "Variable ", vn , " contains only integers but is not type 'integer'. If you intend it as an index variable, you should as.integer() it before passing to ulam." ) )
                        }
                    }# not matrix
                }# all integer
            }# not character
        }
        # pre-scan for any scale() processed variables that might need to be as.numeric()
        for ( i in 1:length(data) ) {
            if ( !is.null( attr(data[[i]],"scaled:center") ) ) {
                # strip the scaling junk
                data[[i]] <- as.numeric( data[[i]] )
            }
        }
    }

    ########################################
    # empty Stan code
    m_funcs <- "" # functions
    m_data <- "" # data
    m_tdata1 <- ""
    m_tdata2 <- ""
    m_pars <- "" # parameters
    m_tpars1 <- "" # transformed parameters, declarations
    m_tpars2 <- "" # transformed parameters, transformation code
    m_model_declare <- "" # local declarations for linear models
    m_model_txt <- "" # general order-trusting code
    m_gq1 <- "" # generated quantities, declarations
    m_gq2 <- "" # gq body

    m_f_reduce_declare <- "" # reduce_sum function
    m_f_reduce <- ""
    if ( threads > 1 & cmdstan==FALSE ) stop( "threads > 1 currently requires cmdstan=TRUE" )
    if ( threads > 1 & log_lik==TRUE ) stop( "threads > 1 currently not compatible with log_lik=TRUE." )

    # functions for manipulating code blocks
    model_declare_append <- function( the_text ) {
        if ( threads==1 ) {
            newcode <- concat( m_model_declare , the_text )
            assign( "m_model_declare" , newcode , inherits=TRUE )
        } else {
            # add instead to m_f_reduce, will go into functions
            newcode <- concat( m_f_reduce_declare , the_text )
            assign( "m_f_reduce_declare" , newcode , inherits=TRUE )
        }
    }
    model_body_append <- function( the_text , override=FALSE ) {
        if ( threads==1 | override==TRUE ) {
            newcode <- concat( m_model_txt , the_text )
            assign( "m_model_txt" , newcode , inherits=TRUE )
        } else {
            # add instead to m_f_reduce, will go into functions
            newcode <- concat( m_f_reduce , the_text )
            assign( "m_f_reduce" , newcode , inherits=TRUE )
        }
    }



    indent <- "    " # 4 spaces
    inden <- function(x) paste( rep(indent,times=x) , collapse="" )

    if ( missing(pars_omit) ) pars_omit <- c()

    flist.orig <- flist
    #flist <- flist_untag(flist) # converts <- to ~ and evals to formulas

    if ( missing(start) ) start <- list()

    # inverse link list
    inverse_links <- list(
        exp = 'log',
        log = 'exp',
        logit = 'inv_logit',
        logodds = 'inv_logit',
        cloglog = 'inv_cloglog'
    )

    # stan variable types and their dimensions
    # these can be used in explicit type declarations preceding symbols, followed by :
    # e.g. vector[10]:z
    stan_types <- list(
        real = 1,
        vector = 1,
        matrix = 2
    )

    declare_templates <- list(
        real = list( patt="real<CONSTRAINT> NAME[DIM1]" , dims=1 ),
        vector = list( patt="vector<CONSTRAINT>[DIM1] NAME[DIM2]" , dims=2 ),
        row_vector = list( patt="row_vector<CONSTRAINT>[DIM1] NAME[DIM2]" , dims=2 ),
        matrix = list( patt="matrix<CONSTRAINT>[DIM1,DIM2] NAME[DIM3]" , dims=3 ),
        int = list( patt="int<CONSTRAINT> NAME[DIM1]" , dims=1 ),
        int_array = list( patt="int<CONSTRAINT> NAME[DIM1,DIM2]" , dims=2 ),
        corr_matrix = list( patt="corr_matrix<CONSTRAINT>[DIM1] NAME[DIM2]" , dims=2 ),
        cholesky_factor_corr = list( patt="cholesky_factor_corr<CONSTRAINT>[DIM1] NAME[DIM2]" , dims=2 ),
        ordered = list( patt="ordered<CONSTRAINT>[DIM1] NAME[DIM2]" , dims=2 ),
        simplex = list( patt="simplex[DIM1] NAME[DIM2]" , dims=2 )
    )

    # binary operators that may appear in linear models
    operations <- c( "~" , "+" , "-" , "*" , "/" , "^" , ":" , "|" , "<-" , "[" )

    # build name lists for distributons, so we can exclude these as symbols later
    dist_names_stan <- unique( sapply( 1:length(distribution_library) , function(L) distribution_library[[L]]$Stan_name ) )
    dist_names_R <- unique( sapply( 1:length(distribution_library) , function(L) distribution_library[[L]]$R_name ) )

    symbols <- list() # data/parameters/variables
    indices <- list() # 'N' data variables to insert
    link_funcs <- list() # link functions built from linear models
    impute_bank <- list() # info for imputation of continuous missing values
    marginalize_bank <- list() # info for marginalization of discrete missing values 
    macros_processed <- c()

    #################################################
    # worker functions

    get_all_symbols <- function( X , strip=TRUE ) {
        # used to get vector of names of symbols in expression X
        # typically used to parse linear models
        # but can parse entire formula lines, stripping links, operators, type declarations
        symbols <- c()
        if ( class(X)=="name" ) {
            # only a symbol here
            X <- as.list(X)
        }
        for ( i in 1:length(X) ) {
            if ( class( X[[i]] )=="call" || class( X[[i]] )=="(" ) {
                # function call, so use recursion to parse it
                nested_symbols <- get_all_symbols( X[[i]] , strip=strip )
                symbols <- c( symbols , nested_symbols )
            }
            if ( class( X[[i]] )=="name" ) {
                if ( strip==TRUE ) {
                    if ( !(as.character(X[[i]]) %in% c( operations , names(inverse_links) , names(stan_types) , dist_names_stan , dist_names_R , "c" ) ) ) {
                        # a symbol, not a binary operator
                        symbols <- c( symbols , as.character(X[[i]]) )
                    }
                } else {
                    # no filtering, just strip out operations etc
                    if ( !(as.character(X[[i]]) %in% c( operations , names(stan_types) , "c" ) ) )
                        symbols <- c( symbols , as.character(X[[i]]) )
                }
            }
        }#i
        return( unique(symbols) )
    }

    detect_left_type <- function( s , ff ) {
        if ( as.character( ff[[1]] ) %in% c('<-','<<-') ) {
            # must be a local variable
            return( "local" )
        } else {
            # now check left side
            if ( any(s %in% names(data)) ) {
                # observed variable
                return( "data" )
            } else {
                return( "par" )
            }
        }
        # if got this far, something wrong
        return( "error" )
    }

    get_left_symbol <- function( f ) {
        # extracts left symbol name, handling explicit dims with : and conditioning with | and indexing with []
        symbol <- NULL
        if ( class( f[[2]] ) == "call" ) {
            # some function here instead of a raw name
            symbol <- as.character( f[[2]] ) # default processing, works for link functions
            # special processing follows
            if ( as.character( f[[2]][[1]] )==':' ) {
                symbol <- as.character( f[[2]][[3]] )
                dims <- as.character(f[[2]][[2]] )
                attr(symbol,"dims") <- as.list( dims[2:length(dims)] )
                if ( symbol[1]=="[" ) {
                    # also has a grouping variable (in array position)
                    sa <- list( dims[2] , dims[3] , symbol[3] )
                    symbol <- symbol[-c(1,3)]
                    attr(symbol,"dims") <- sa
                }
                if ( symbol[1]=="[[" ) {
                    # index for specific assignment
                    # just need to prune the symbol name down
                    sa <- attr(symbol,"dims")
                    symbol <- symbol[-c(1,3)]
                    attr(symbol,"dims") <- sa
                }
            }
            if ( as.character( f[[2]][[1]] )=='|' ) {
                symbol <- as.character( f[[2]][[2]] )
                attr(symbol,"iff") <- f[[2]][[3]]
            }
            if ( as.character( f[[2]][[1]] )=='[' ) {
                # should be a length or symbol to figure length out
                symbol <- as.character( f[[2]][[2]] )
                attr(symbol,"dims") <- f[[2]][[3]]
            }
            if ( as.character( f[[2]][[1]] )=='[[' ) {
                # a single element of a vector, used to assign different priors to each element
                # can ignore the [[n]] until distribution is composed (in compose_distribution)
                symbol <- as.character( f[[2]][[2]] )
            }
            if ( as.character( f[[2]][[1]] ) == 'c' ) {
                # vector construction like c(a,b)
                # should parse fine into list like "c" "a" "b"
                # composition function later handles multiple implied declarations
                symbol <- as.character( f[[2]] )
                symbol[1] <- "VECTOR"
            }
            if ( as.character( f[[2]][[1]] ) == '>' ) {
                # block tag
                # 'save>' places local assignment in transformed parameters
                # 'gq>' places in generated quantities
                the_block <- deparse( f[[2]][[2]] )

                # need to send it down recursively to process other features
                ff <- f
                ff[[2]] <- ff[[2]][[3]]
                symbol <- get_left_symbol( ff )
                attr(symbol,"block") <- the_block
            }
        } else {
            # a raw name
            symbol <- as.character( f[[2]] )
        }
        # check for link functions
        if ( length(symbol) == 2 ) {
            if ( symbol[1] %in% names(inverse_links) ) {
                new_symbol <- symbol[2]
                mostattributes(new_symbol) <- attributes(symbol)
                attr(new_symbol,"link") <- symbol[1]
                symbol <- new_symbol
            }
        }
        if ( symbol[1] == 'c' && length(symbol)>1 ) {
            # vector construction like c(a,b)
            new_symbol <- symbol
            new_symbol[1] <- "VECTOR"
            mostattributes(new_symbol) <- attributes(symbol)
            symbol <- new_symbol
        }
        # result
        return( symbol )
    }

    get_dist_template <- function( the_dist ) {
        found <- FALSE
        for ( j in 1:length(distribution_library) ) {
            if ( distribution_library[[j]]$Stan_name==the_dist || distribution_library[[j]]$R_name==the_dist || names(distribution_library)[j]==the_dist ) {
                template <- distribution_library[[ j ]]
                found <- TRUE
                break
            }
        }
        if ( found==FALSE ) stop( concat( "No template for distribution '",the_dist,"'" ) )
        return( template )
    }

    # gsub( "NAME" , "par" , "int NAME" , fixed=TRUE )
    compose_declaration <- function( symbol_name, symbol_list , reduce=NULL ) {
        # the_dims should be a list like: [[1]] "vector" [[2]] 5

        if ( class( symbol_list$dims ) == "character" )
            # just a stan type, so assume scalar
            symbol_list$dims <- list( symbol_list$dims , 1 )

        if ( class( symbol_list$dims ) == "name" ) {
            # has a grouping variable, so need to fetch its UNIQUE length
            ul <- length( unique( data[[ as.character(symbol_list$dims) ]] ) )
            symbol_list$dims <- list( "vector" , ul )
        }

        template <- declare_templates[[ symbol_list$dims[[1]] ]]
        if ( is.null(template) ) stop( concat( "Declaration template not found: ", symbol_list$dims[[1]] ) )
        n_dims <- template$dims
        out <- template$patt
        out <- gsub( "NAME" , symbol_name , out , fixed=TRUE )
        constraint_txt <- ""
        if ( !is.null(symbol_list$constraint) ) {
            if ( !is.na(symbol_list$constraint) ) {
                constraint_txt <- concat( "<" , symbol_list$constraint , ">" )
            }
        }
        out <- gsub( "<CONSTRAINT>" , constraint_txt , out , fixed=TRUE )

        for ( i in 1:n_dims ) {
            Z <- concat( "DIM" , i )
            if ( (i+1) <= length( symbol_list$dims ) ) {
                d <- symbol_list$dims[[i+1]]
                if ( !is.null( data[[ as.character(d) ]] ) ) {
                    # either a grouping variable OR an index length (scalar)
                    if ( length( data[[ as.character(d) ]] ) > 1 ) {
                        # a grouping variable, so count number of unique values
                        d <- length(unique( data[[ as.character(d) ]] ))
                    }
                    # if an index, don't need to do anything, converted to text later
                }
                dd <- as.character( d )
                if ( !is.null(reduce) ) dd <- concat("size(",reduce,")")
                # don't add dims if scalar and in array dim spot (at end)
                if ( i==n_dims && d==1 ) {
                    dd <- ""
                    Z <- concat("[",Z,"]")
                }
                out <- gsub( Z , dd , out , fixed=TRUE )
            } else {
                # erase optional array dim
                # could be either [DIM2] or [DIM1,DIM2]
                out <- gsub( concat("[",Z,"]") , "" , out , fixed=TRUE )
                out <- gsub( concat(",",Z,"]") , "]" , out , fixed=TRUE )
            }
        }
        return( out )
    }

    compose_distibution <- function( left_symbol , fline , as_log_lik=FALSE ) {
        # this function builds text line to insert into model block
        # form: symbol ~ distibution( par1 , par2 , ... )
        # typically needs to no modification
        # any modification done by template, allowing some crazy macro stuff

        # is there a [[n]] on the left side?
        if ( class( fline[[2]] )=="call" ) {
            if ( as.character( fline[[2]][[1]] )=="[[" ) {
                # augment left symbol
                the_index <- as.character( fline[[2]][[3]] )
                left_symbol <- concat( left_symbol , "[" , the_index , "]" )
            }
            # check also after an explicit dim ':'
            if ( as.character( fline[[2]][[1]] )==":" ) {
                if ( class( fline[[2]][[3]] )=="call" ) {
                    if ( as.character( fline[[2]][[3]][[1]] )=="[[" ) {
                        # augment left symbol
                        the_index <- as.character( fline[[2]][[3]][[3]] )
                        left_symbol <- concat( left_symbol , "[" , the_index , "]" )
                    }
                }
            }
        }

        right <- as.character( fline[[3]] ) # will be a character vector including pars
        left <- left_symbol # any bracket should have been stripped away already
        template <- get_dist_template( right[1] )
        if ( is.null( template$build ) )
            build_func <- distribution_library[['DEFAULT']]$build
        else
            build_func <- template$build
        environment( build_func ) <- parent.env(environment())
        if ( is.null( template$build ) )
            out <- do.call( build_func , list( left , fline[[3]] , template$Stan_name , template$Stan_suffix , as_log_lik=as_log_lik , vectorized=template$vectorized ) , quote=TRUE )
        else
            out <- do.call( build_func , list( left , fline[[3]] , as_log_lik=as_log_lik ) , quote=TRUE )
        # done
        return( out )
    }

    # insert [i] on each data symbol
    # this is a recursive function, drilling down in expression tree
    # need to be careful not to double-up indexes like theta[i][group[i]]
    nest_indexes <- function( f , data_symbols ) {
        
        for ( i in 1:length(f) ) {
            if ( class( f[[i]] )=="name" ) {
                # could be a data (or local) var
                if ( as.character(f[[i]]) %in% data_symbols ) {
                    # check if this already has a `[` next to it
                    if ( i==2 )
                        if ( as.character(f[[1]])=="[" )
                            next # move on to next symbol, bc this already has [ after it
                    # check if length 1, so a constant that doesn't need [i]
                    if ( length( data[[ as.character(f[[i]]) ]] )==1 ) next
                    # convert to a nested call
                    f[[i]] <- call( "[" , f[[i]] , quote(i) )
                }
            } else {
                if ( class( f[[i]] )=="call" || class( f[[i]] )=="(" ) {
                    # nested structure
                    # need to drill down
                    f[[i]] <- nest_indexes( f[[i]] , data_symbols )
                }
            }
        }#i

        return( f )
    }

    # recursively search formula for '!Y' and replace with new symbol
    # used by custom template
    hunt_bangY <- function( f , r ) {
        if ( class( f )=="call" || class( f )=="(" ) {
            if ( as.character( f[[1]] )=="!" ) {
                # bang found - now check argument
                arg <- deparse(f[[2]])
                if ( arg=="Y" ) {
                    # got it - now replace entire f (replace ! too)
                    f <- as.name( r )
                }
            } else {
                # need to drill down
                for ( i in 1:length(f) ) f[[i]] <- hunt_bangY( f[[i]] , r )
            }
        }
        return( f )
    }

    compose_assignment <- function( symbol , fline , reduce=NULL ) {
        # compose local assignment (usually a linear model) for model block
        # need to observe need for loop over assignment and insert [i] after the data symbols on right-hand side
        # when assignment is <<- instead then no loop (vectorized assignment)
        # need to check also for some operator conversions R -> Stan:
        #     %*% -> *
        # can do this in text as last step? not elegant.
        the_dims <- symbols[[symbol]]$dims
        N <- NA
        if ( class(the_dims)=="character" ) {
            # just a type label
            N <- 1
        } else {
            # should be a list
            N <- the_dims[[2]]
            if ( the_dims[[1]]=="matrix" ) {
                # matrix local, so no loop? have to guess intent is matrix multiplication on right
                N <- NA
            }
        }
        # build loop text, if any
        loop_txt <- ""
        local_indent <- indent
        symbol_txt <- symbol
        rline <- fline[[3]]
        rline_link <- rline
        assignment_symbol <- as.character( fline[[1]] )
        if ( assignment_symbol == "<<-" ) N <- NA # mark to vectorize

        # link version needs [TRUE,] nesting done before [i], otherwise end up with [TRUE,i] on data variables
        rline_link <- nest_commas( rline_link )

        if ( !is.na(N) ) {
            if ( N > 1 ) {
                # if in reduce context, need size(outcome) as length
                # not likely to know name of outcome var at this point, so insert XREDUCEOUTX placeholder - will replace later with known outcome
                N_insert <- N
                if ( reduce==TRUE ) N_insert <- concat("size(XREDUCEOUTX)")
                loop_txt <- concat( indent , "for ( i in 1:" , N_insert , " ) {\n" )
                local_indent <- inden(2)
                # insert [i] to each data variable on right-hand side
                # do this by nesting call in expression tree, then can use deparse to convert all to text
                right_symbols <- get_all_symbols( rline )
                data_symbols <- right_symbols[ which( right_symbols %in% names(data) ) ]
                for ( k in 1:length(symbols) ) { 
                    # also add any locals that are lower in formula list
                    # but if local is matrix type, then DO NOT add [i], so leave out of this list
                    if ( symbols[[k]]$type=="local" ) {
                        if ( symbols[[k]]$dims[[1]] != "matrix" )
                            data_symbols <- c( data_symbols , names(symbols)[k] )
                    }
                }

                rline <- nest_indexes( rline , data_symbols )
                rline_link <- nest_indexes( rline_link , data_symbols )
                # add [i] to left symbol
                symbol_txt <- concat( symbol_txt , "[i]" )
            }
        }

        # save parsed version to build link function later
        # symbols[[symbol]]$right_parsed <- rline
        new_symbols <- symbols
        new_symbols[[symbol]]$right_iparsed <- rline_link
        new_symbols[[symbol]]$link_N <- ifelse( is.na(N) , 1 , N )
        assign( "symbols" , new_symbols , inherits=TRUE )

        # compose text - should have [i] in it already
        # check for transpose t() function and convert to Stan's trailing '
        transpose_flag <- FALSE
        if ( class(rline)=="call" ) {
            if ( as.character( rline[[1]] )=="t" ) {
                # transpose found
                transpose_flag <- TRUE
                # cut it off and add it on end after deparse
                rline[[1]] <- as.name("(") # replace function name
            }
        }
        
        right_txt <- deparse( rline , width.cutoff = 500L )
        if ( length(right_txt)>1 ) {
            # lines got broken! so concatenate with line breaks
            right_text <- paste( right_txt , collapse="\n" )
        }
        # if reduce context, need [start+i-1] in place of [i]
        if ( reduce==TRUE )
            right_txt <- gsub( "[i]" , "[start+i-1]" , right_txt , fixed=TRUE )
        if ( transpose_flag==TRUE ) {
            right_txt <- concat( right_txt , "'" )
        }

        out <- concat( local_indent , symbol_txt , " = " , right_txt , ";\n" )
        if ( !is.na( symbols[[symbol]]$link ) ) {
            inv_link <- inverse_links[[ symbols[[symbol]]$link ]]
            if ( is.null( inv_link ) ) stop( concat( "No inverse link found for link '" , symbols[[symbol]]$link , "'" ) )
            out <- concat( out , local_indent , symbol_txt , " = " , inv_link , "(" , symbol_txt , ");\n" )
        }
        if ( !is.na(N) )
            if ( N > 1 ) {
                # prepend and close loop
                out <- concat( loop_txt , out , indent , "}\n" )
            }

        # finally check for any function conversions R -> Stan
        R_to_Stan_f <- list( 
            " %*% " = " * " ,
            " as.vector(" = " to_vector("
        )
        for ( j in 1:length(R_to_Stan_f) )
            out <- gsub( names(R_to_Stan_f)[j] , R_to_Stan_f[[j]] , out , fixed=TRUE )

        # all done
        return( out )
    }

    register_data_var <- function( var_name ) {

        # get dims from data list
        the_dims <- dim( data[[var_name]] )
        if ( is.null(the_dims) ) the_dims <- length( data[[var_name]] )
        
        # try to determine Stan type from class
        stan_type <- "real"
        if ( class( data[[var_name]] )[1]=="integer" ) stan_type <- "int"
        if ( inherits( data[[var_name]] , "matrix" ) ) {
            stan_type <- "matrix"
            the_dims <- list( stan_type , the_dims[1] , the_dims[2] )
        } else {
            # numeric or integer
            the_dims <- list( stan_type , the_dims )
        }
        # check for vector and call it a vector
        if ( stan_type=="real" ) {
            if ( the_dims[[2]] > 1 ) the_dims[[1]] <- "vector"
        }
        # check for integer array
        if ( class( data[[var_name]] )[1]=="array" ) {
            if ( class( data[[var_name]][1] )[1]=="integer" ) {
                the_dims <- list( "int_array" , dim( data[[var_name]] ) )
            } else {
                the_dims <- list( "real" , dim( data[[var_name]] ) )
            }
        }
        # return list suitable to insert in symbols list
        return( list( name=var_name , type="data" , dims=the_dims , constraint=NA ) )
    }

    # adds those [TRUE,i] things in link functions
    nest_commas <- function( f ) {

        if ( class( f )=="call" || class( f )=="(" ) {
            if ( as.character( f[[1]] )=="[" ) {
                # bracket call, so insert comma (TRUE)
                the_call <- as.list( f )
                f <- as.call( unlist( list( the_call[1:2] , quote(TRUE) , the_call[3:length(the_call)] ) , recursive=FALSE ) )
            } else {
                # need to drill down
                for ( i in 1:length(f) )
                    f[[i]] <- nest_commas( f[[i]] )
            }
        }

        return( f )
    }

    # search expression tree recursively and replace pattern with x
    symbol_gsub <- function( f , pattern , x ) {
        for ( i in 1:length(f) ) {
            if ( class( f[[i]] )=="name" ) {
                # could be pattern
                if ( as.character(f[[i]]) == pattern ) {
                    f[[i]] <- as.name( x )
                }
            } else {
                if ( class( f[[i]] )=="call" || class( f[[i]] )=="(" ) {
                    # nested structure - need to drill down
                    f[[i]] <- symbol_gsub( f[[i]] , pattern , x )
                }
            }
        }#i
        return( f )
    }

    compose_link_func <- function( right , the_link , N ) {
        # need to insert leading commas in bracket expressions like [group] -> [TRUE,group]
        # this allows link function to process over all samples (samples in first index position)
        # special problem: some local variables are functions of parameters, and so must be treated like parameters when embedded in other local variables
        # will try to solve that inside link() at run time
        #right_prime <- nest_commas( right )
        right_prime <- right
        the_inv_link <- inverse_links[[ the_link ]]
        out <- list( func=right_prime , link=the_link , inv_link=the_inv_link , N=N )
        return(out)
    }

    process_macro <- function( the_macro , flist , n ) {
        # call build function
        build_func <- macro_library[[ the_macro ]]$build
        if ( !is.null( build_func ) ) {
            environment( build_func ) <- parent.env(environment())
            flist <- do.call( build_func , list( flist , n ) , quote=TRUE )
        }
        return( flist )
    }

    #################################################
    # pre-parse scan for NA values and process proper substitutions or errors
    lf <- length(flist)
    for ( i in lf:1 ) {
        all_symbols <- get_all_symbols( flist[[i]] )
        idx <- which( all_symbols %in% names(data) )
        if ( length(idx)>0 ) {
            # found some data symbols
            # check for NA values
            for ( j in idx ) {
                if ( any( is.na( data[[ all_symbols[j] ]] ) ) ) {
                    var <- all_symbols[j]
                    n_miss <- sum( is.na( data[[ var ]] ) )
                    # check whether already in impute bank
                    if ( is.null( impute_bank[[ var ]] ) ) {
                        # on left side? if so, then check the type and set up imputation
                        the_left <- get_left_symbol( flist[[i]] )
                        if ( the_left[1]==var ) {
                            the_type <- class( data[[ var ]] )

                            if ( the_type=="numeric" ) {
                                # continuous --- set up merge variable
                                var_merge <- concat( var , "_merge" ) # for linear models
                                var_impute <- concat( var , "_impute" ) # pars for missing values
                                # replace left side with var_merge symbol
                                flist[[i]][[2]] <- parse( text=var_merge )[[1]]
                                # insert merge_missing macro below this line in the formula list
                                merge_line <- concat( var_merge , " <- merge_missing( " , var , " , " , var_impute , " )" )
                                merge_line <- parse( text=merge_line )[[1]]
                                flist <- append( flist , merge_line , i )
                                # add impute_bank entry
                                impute_bank[[ var ]] <- list( type="real" , n=n_miss , merge=var_merge )
                                # message
                                if ( messages==TRUE )
                                    message( concat( "Found " , n_miss , " NA values in " , var , " and attempting imputation." ) )
                            }

                            if ( the_type=="integer" ) {
                                # discrete --- set up marginalization mixture
                            }

                        }# on left
                        else {
                            # on right - but not already in impute bank
                            # could be a merge_missing call, so check
                            right_token <- as.character( flist[[i]][[3]] )[1]
                            if ( right_token != "merge_missing" )
                                if ( messages==TRUE )
                                    message( concat( "Found " , n_miss , " NA values in " , var , ". Not sure what to do with these, and they might precent the model from running." ) )
                        }

                    } else {
                        # already in impute_bank
                        # check if we are on right-hand side and need to replace with merge symbol
                        right_symbols <- get_all_symbols( flist[[i]][[3]] )
                        if ( var %in% right_symbols ) {
                            # replace with merge symbol
                            var_merge <- impute_bank[[ var ]]$merge
                            flist[[i]][[3]] <- symbol_gsub( flist[[i]][[3]] , var , var_merge )
                        }# on right
                    }
                }
            }#j
        }
    }#i

    #################################################
    # check for any macros
    # macros are symbols in the macro_library argument (by default ulam_macros from rethinking namespace)
    # once a symbol is recognized, its build() function is called, passing calling environment
    # any code in 'functions' slot is added to Stan functions block
    # as a result, one way to add custom distributions is to declare them in both the macros and distributions libraries
    # macro can add the Stan function
    # distribution template handles the arguments and constraints

    # formula list processed from bottom to top
    # macros can add arbitrary Stan code to any block
    # macros can *add* extra formula lines *after* own line without messing things up
    # macros can modify *any* formula line
    lf <- length(flist)
    flist.new <- flist
    for ( i in lf:1 ) {
        all_symbols <- get_all_symbols( flist.new[[i]] , strip=FALSE )
        if ( any( all_symbols %in% names(macro_library) ) ) {
            # hit
            idx <- which( all_symbols %in% names(macro_library) )
            # for each hit, process the macro
            for ( k in 1:length(idx) ) {
                the_macro <- all_symbols[ idx[k] ]
                if ( !(the_macro %in% macros_processed) ) {
                    # check for code to insert in Stan code
                    if ( !is.null( macro_library[[ the_macro ]]$functions ) ) {
                        m_funcs <- concat( m_funcs , "\n" , macro_library[[ the_macro ]]$functions , "\n" )
                    }
                    if ( !is.null( macro_library[[ the_macro ]]$transformed_data_declare ) ) {
                        m_tdata1 <- concat( m_tdata1 , "\n" , macro_library[[ the_macro ]]$transformed_data_declare , "\n" )
                    }
                    if ( !is.null( macro_library[[ the_macro ]]$transformed_data_body ) ) {
                        m_tdata2 <- concat( m_tdata2 , "\n" , macro_library[[ the_macro ]]$transformed_data_body , "\n" )
                    }
                }
                # do the build function, if it exists
                flist.new <- process_macro( the_macro , flist.new , i )
                # add the_macro to macros_processed list so we don't add any code blocks again
                macros_processed <- unique( c( macros_processed , the_macro ) )
            }#k
        }
    }#i
    flist <- flist.new

    #################################################
    # the action loop
    # for each line in formula list, parse symbols and build corresponding code lines
    # some parsing is pending later lines, due to ambiguities
    # goal is to focus on left-hand symbols, as each should have definition line (prior or likelihood or assignment)
    # but some data variables will only appear on right-hand side of distributions/equations, so have to scan right-hand sides for data variables to declare

    # first pass just to ID left-hand symbols and build a graph of connections
    # this will let us decide some things on second pass (which actively parses)

    symbol_graph <- NULL
    symbol_lines <- NULL
    nobs_save <- 0
    # how to read the symbol graph:
    # rows contain the columns
    # columns are contained in the rows

    # find the nodes
    for ( i in 1:length(flist) ) {
        left_symbol <- get_left_symbol( flist[[i]] )
        left_type <- detect_left_type( left_symbol , flist[[i]] )
        if ( length(left_symbol) > 1 ) {
            # a vector of parameters?
            if ( left_symbol[1] == "VECTOR" ) {
                new_left_symbol <- left_symbol
                new_left_symbol <- new_left_symbol[-1]
                mostattributes(new_left_symbol) <- attributes(left_symbol)
                left_symbol <- new_left_symbol
            }
        }
        old_dim <- nrow(symbol_graph)
        if ( i==1 ) old_dim <- 0
        new_dim <- old_dim + length(left_symbol)
        new_graph <- matrix( 0 , nrow=new_dim , ncol=new_dim )

        # manage margin names
        if ( i > 1 ) {
            new_graph[ 1:old_dim , 1:old_dim ] <- symbol_graph
            colnames(new_graph) <- c( colnames(symbol_graph) , rep(NA,length(left_symbol)) )
            rownames(new_graph) <- colnames(new_graph)
        } else {
            # init names
            rownames( new_graph ) <- rep("XX",new_dim)
            colnames( new_graph ) <- rep("XX",new_dim)
        }
        for ( j in 1:length(left_symbol) ) {
            rownames( new_graph )[old_dim+j] <- left_symbol[j]
            colnames( new_graph )[old_dim+j] <- left_symbol[j]
        }
        symbol_graph <- new_graph

        # store formula line number --- different than graph margin index when more than one left symbol on a line
        for ( j in 1:length(left_symbol) ) {
            symbol_lines[[ left_symbol[j] ]] <- i
        }
    }#i
    # scan for connections now
    for ( i in 1:nrow(symbol_graph) ) {
        left_symbol <- rownames(symbol_graph)[i]
        # scan symbols on right-hand side
        right_symbols <- NULL
        if ( i <= length(symbol_lines) )
            right_symbols <- get_all_symbols( flist[[ symbol_lines[[i]] ]][[3]] )
        if ( length(right_symbols)>0 ) {
            for ( j in 1:length(right_symbols) ) {
                if ( right_symbols[j] %in% colnames(new_graph) ) {
                    # flag all left symbols
                    idx_left <- which( colnames(new_graph) %in% left_symbol )
                    idx_right <- which( colnames(new_graph)==right_symbols[j] )
                    symbol_graph[ idx_left , idx_right ] <- 1
                }
            }#j
        }
    }#i

    thflag <- FALSE
    reduce_outcome <- NULL

    # second pass - build code
    # this pass goes backwards
    for ( i in length(flist):1 ) {

        left_symbol <- get_left_symbol( flist[[i]] )
        left_type <- detect_left_type( left_symbol , flist[[i]] )

        if ( length(left_symbol) > 1 ) {
            # a vector of parameters?
            if ( left_symbol[1] == "VECTOR" ) {
                new_left_symbol <- left_symbol
                new_left_symbol <- new_left_symbol[-1]
                mostattributes(new_left_symbol) <- attributes(left_symbol)
                left_symbol <- new_left_symbol
            }
        }

        if ( left_type=="par" ) {

            # need to find dims
            the_dims <- NA
            constraint <- NA
            the_dist <- as.character( flist[[i]][[3]] ) # will be a character vector including pars
            template <- get_dist_template( the_dist[1] )

            constraint <- template$constraints[1] # implied constraint
            # check for custom contraint
            if ( !missing(constraints) ) {
                for ( kk in 1:length(constraints) ) {
                    if ( names(constraints)[kk]==left_symbol[1] ) {
                        constraint <- constraints[[kk]] # should be text
                    }
                }#kk
            }
            
            if ( !is.null( attr(left_symbol,"dims") ) ) 
                # dims must have been explicit and already parsed by get_left_symbol
                # also check for any implied constraints
                the_dims <- attr(left_symbol,"dims")
            else {
                # try getting implied dims from right side distribution
                # get default dim of outcome to start
                the_dims <- template$dims[1]
            }

            # check for data vars on right side and check types as well
            # these are registered further down, but set type now when have template info
            right <- the_dist[-1]
            for ( k in 1:length(right) ) {
                if ( right[k] %in% names(data) ) {
                    # data so check type
                    need_type <- template$dims[k+1]
                    if ( need_type %in% c("real","vector") ) {
                        if ( class(data[[ right[k] ]])[1] != "numeric" )
                            data[[ right[k] ]] <- as.numeric(data[[ right[k] ]])
                    }
                }
            }#k

            # store info
            # might be multiple parameters in a vector
            for ( j in 1:length(left_symbol) ) {
                # check whether this symbol already declared as local in previous line
                # if so, don't mess up its type and dims here!
                do_add <- TRUE
                if ( ( left_symbol[j] %in% names(symbols) ) ) {
                    # check if type=='local' or previous dimming was explicit (: operator)
                    if ( symbols[[ left_symbol[j] ]]$type=="local" ) do_add <- FALSE
                }
                if ( do_add==TRUE )
                    symbols[[ left_symbol[j] ]] <- list( type=left_type , dims=the_dims , constraint=constraint )
            }#j

            # build text for model block
            built <- compose_distibution( left_symbol , flist[[i]] )
            #m_model_txt <- concat( m_model_txt , built )
            model_body_append( built , TRUE )

            # apply constraints on any parameters on right-hand side
            # would be better to do this inside compose function, but outside avoids scope issues
            # also a good time to check/update dims on these variables
            the_dist <- as.character( flist[[i]][[3]] )
            template <- get_dist_template( the_dist[1] )
            for ( j in 2:length(flist[[i]][[3]]) ) {
                # 1st index should just be distribution name
                # following positions are parameters, possibly complex expressions
                # when find a symbol, check if it is already in symbol list
                # if in list, then apply implied constraint from distribution template
                if ( class( flist[[i]][[3]][[j]] )=="name" ) {
                    par_name <- as.character( flist[[i]][[3]][[j]] )
                    # only overwrite constraints if no meaningful constraints there already
                    if ( par_name %in% names(symbols) ) {
                        # constraint?
                        if ( is.null(symbols[[par_name]]$constraint) )
                            symbols[[par_name]]$constraint <- template$constraints[j]
                        else
                            if ( is.na(symbols[[par_name]]$constraint) )
                                symbols[[par_name]]$constraint <- template$constraints[j]
                        # dimensions?
                    }
                }
            }#j

        }#par

        if ( left_type=="data" & sample_prior==FALSE ) {

            # first check for implied Stan type for left var
            the_dist <- as.character( flist[[i]][[3]] ) # will be a character vector including pars
            template <- get_dist_template( the_dist[1] )
            out_type <- template$dims[1] # Stan type

            # now process each var on left side
            for ( j in 1:length(left_symbol) ) {
                if ( !is.na(out_type) ) {
                    if ( out_type=="real" ) 
                        data[[ left_symbol[j] ]] <- as.numeric(data[[ left_symbol[j] ]])
                    if ( out_type=="int" ) 
                        data[[ left_symbol[j] ]] <- as.integer(data[[ left_symbol[j] ]])
                }
                symbols[[ left_symbol[j] ]] <- register_data_var( left_symbol[j] )
            }
            # check for data vars on right side and check types as well
            # these are registered further down, but set type now when have template info
            right <- the_dist[-1]
            for ( k in 1:length(right) ) {
                if ( right[k] %in% names(data) ) {
                    # data so check type
                    need_type <- template$dims[k+1]
                    if ( need_type %in% c("real","vector") ) {
                        if ( class(data[[ right[k] ]]) != "numeric" )
                            data[[ right[k] ]] <- as.numeric(data[[ right[k] ]])
                    }
                }
            }#k

            # build text for model block
            # destined for reduce_sum? if so format right
            if ( threads > 1 ) {
                thflag <- "reduce"
                reduce_outcome <- left_symbol
            }
            built <- compose_distibution( left_symbol , flist[[i]] , as_log_lik=thflag )
            #m_model_txt <- concat( m_model_txt , built )
            model_body_append( built )
            # add log_lik to generated quantities?
            # currently just uses first outcome --- for multiple outcome models, will need new code
            if ( i==1 && log_lik==TRUE ) {
                # add by default
                built <- compose_distibution( left_symbol , flist[[i]] , as_log_lik=TRUE )
                m_gq2 <- concat( m_gq2 , built )
                N <- symbols[[left_symbol]]$dims[[2]]
                m_gq1 <- concat( m_gq1 , indent , "vector[" , N , "] log_lik;\n" )
                # save N to attr so nobs/compare can get it later
                nobs_save <- N
            }
        }#data

        if ( left_type=="local" ) {
            # local model block variable
            # dims come from data symbol length on right side or from later distribution that uses this local
            the_dims <- NA
            if ( !is.null( attr(left_symbol,"dims") ) ) 
                the_dims <- attr(left_symbol,"dims")
            else {
                # search right side, find length of a data variable
                right_symbols <- get_all_symbols( flist[[i]][[3]] )
                if ( length(right_symbols)==0 ) stop( concat( "Tried to parse a linear model with ZERO symbols:", deparse(flist[[i]]) ) )
                hits <- which( right_symbols %in% names( data ) )
                if ( length(hits) > 0 ) {
                    first_hit <- right_symbols[ hits[1] ]
                    the_dims <- list( "vector" , length( data[[ first_hit ]] ) )
                } else {
                    # no data var hits
                    # so try to get length from parameters
                    for ( ii in 1:length(right_symbols) ) {
                        if ( right_symbols[ii] %in% names(symbols) ) {
                            if ( symbols[[ii]]$type=="par" ) {
                                # check its dims
                                j <- which( names(symbols)==right_symbols[ii] )
                                if ( any(!is.na( symbols[[j]]$dims )) ) {
                                    zzdims <- 0
                                    if ( length(symbols[[j]]$dims)==1 ) zzdims <- 1
                                    else zzdims <- symbols[[j]]$dims[[2]]
                                    the_dims <- max( the_dims , zzdims , na.rm=TRUE )
                                }
                            }
                        }
                    }#ii
                    # if still NA at this point, then maybe no prior for symbol
                    if ( is.na(the_dims) ) stop(concat("Unable to determine type and dimensions for ",left_symbol[1])," - Do all parameters have priors?")
                    if ( the_dims==1 )
                        the_dims <- "real"
                    else
                        the_dims <- list( "vector" , the_dims )
                }
            }

            # register symbol
            the_link <- NA
            if ( !is.null( attr(left_symbol,"link") ) ) the_link <- attr(left_symbol,"link")
            the_block <- "model"
            if ( !is.null( attr(left_symbol,"block") ) ) the_block <- attr(left_symbol,"block")
            symbols[[ left_symbol ]] <- list( type=left_type , dims=the_dims , link=the_link , block=the_block )

            # build the link function that can be called after sampling
            #link_funcs[[ left_symbol ]] <- compose_link_func( flist[[i]][[3]] , the_link )

        }#local

        # now scan right side for data variables, in case they don't have own likelihoods later
        right_symbols <- get_all_symbols( flist[[i]][[3]] )
        if ( any( right_symbols %in% names(data) ) ) {
            # at least one data variable
            idx <- which( right_symbols %in% names(data) )
            for ( j in idx ) {
                # register each, as long as it hasn't already been registered
                var <- right_symbols[j]
                if ( !( var %in% names(symbols) ) )
                    symbols[[var]] <- register_data_var( var )
            }
        }
        if ( left_type=="local" ) {
            # add assignment to model block
            # do this at end here, because we parsed right-hand side data vars just above
            local_built <- compose_assignment( left_symbol , flist[[i]] , reduce=(threads>1) )
            
            # check for explicit block tag
            the_block <- "model"
            if ( !is.null( attr(left_symbol,"block") ) ) {
                the_block <- attr(left_symbol,"block")
                if ( the_block=="transpars" ) {
                    # transformed parameters
                    m_tpars2 <- concat( m_tpars2 , local_built )
                }
                if ( the_block=="transdata" ) {
                    # transformed data
                    m_tdata2 <- concat( m_tdata2 , local_built )
                }
                if ( the_block=="save" ) {
                    # generated quantities AND model block
                    m_gq2 <- concat( m_gq2 , local_built )
                    #m_model_txt <- concat( m_model_txt , local_built )
                    model_body_append( local_built )
                }
                if ( the_block=="gq" ) {
                    # generated quantities ONLY
                    m_gq2 <- concat( m_gq2 , local_built )
                }
            } else { 
                #m_model_txt <- concat( m_model_txt , local_built )
                model_body_append( local_built )
            }
            
            # add to gq for log_lik, but only when not in transformed pars
            if ( log_lik==TRUE & the_block=="model" ) {
                # also add to generated quantities, so we can compute log_lik, assuming likelihood might contain this local symbol
                m_gq2 <- concat( m_gq2 , local_built )
            }
        }

    }#i

    # build index variables from scanning symbols
    # NOT YET IMPLEMENTED

    # add any remaining data variables named in data list
    # index variables that only appear as left-side indexes e.g.
    if ( declare_all_data==TRUE ) {
        for ( i in names(data) ) {
            if ( !( i %in% names(symbols) ) ) {
                symbols[[ i ]] <- register_data_var( i )
            }
        }
    }

    # add parameter declarations to parameter block
    # add local var declarations to model block
    # add data declarations to data block
    # WE DO THIS BACKWARDS SO PARS/DATA ARE IN ORIGINAL ORDER AS IN FORMULA
    use_pars <- c()
    for ( i in length(symbols):1 ) {
        left_symbol <- names(symbols)[i]
        if ( symbols[[i]]$type=="data" ) {
            Z <- compose_declaration( names(symbols)[i] , symbols[[i]] )
            m_data <- concat( m_data , indent , Z , ";\n" )
            # check for type compatibility on R and Stan sides
            stan_type <- symbols[[i]]$dims
            if ( stan_type[[1]]=="int" )
                data[[ names(symbols)[i] ]] <- as.integer( data[[ names(symbols)[i] ]] )
        }
        if ( symbols[[i]]$type=="par" ) {
            Z <- compose_declaration( names(symbols)[i] , symbols[[i]] )
            m_pars <- concat( m_pars , indent , Z , ";\n" )
        }
        if ( symbols[[i]]$type=="local" ) {

            do_reduce <- NULL
            # use reduce outcome symbol for size when a local declaration in model block - this gives correct size(symbol) in reducer function
            # other blocks need ordinary size declaration
            if ( symbols[[i]]$block=="model" ) do_reduce <- reduce_outcome
            Z <- compose_declaration( names(symbols)[i] , symbols[[i]] , reduce=do_reduce )

            # check for explicit block tag
            the_block <- symbols[[i]]$block
            if ( the_block=="transpars" ) {
                # transformed parameters
                m_tpars1 <- concat( m_tpars1 , indent , Z , ";\n" )
                use_pars <- c( use_pars , left_symbol )
            }
            if ( the_block=="transdata" ) {
                # transformed data
                m_tdata1 <- concat( m_tdata1 , indent , Z , ";\n" )
            }
            if ( the_block=="gq" ) {
                # generated quantities ONLY
                m_gq1 <- concat( m_gq1 , indent , Z , ";\n" )
                use_pars <- c( use_pars , left_symbol )
            }
            if ( the_block=="save" ) {
                # generated quantities AND model block
                m_gq1 <- concat( m_gq1 , indent , Z , ";\n" )
                #m_model_declare <- concat( m_model_declare , indent , Z , ";\n" )
                model_declare_append( concat( indent , Z , ";\n" ) )
                use_pars <- c( use_pars , left_symbol )
            }
            if ( the_block=="model" ) {
                # model block (default)
                #m_model_declare <- concat( m_model_declare , indent , Z , ";\n" )
                model_declare_append( concat( indent , Z , ";\n" ) )
            }

            if ( log_lik==TRUE & the_block %in% c("model","save") ) {
                # also add to generated quantities
                m_gq1 <- concat( m_gq1 , indent , Z , ";\n" )
                pars_omit <- c( pars_omit , left_symbol ) # mark to omit from returned samples
            }
            # compose link function
            if ( the_block=="model" )
            link_funcs[[ left_symbol ]] <- compose_link_func( symbols[[left_symbol]]$right_iparsed , symbols[[left_symbol]]$link , N=symbols[[left_symbol]]$link_N )
        }
    }#i

    # build Stan code
    blockify <- function(x,header,footer) {
        if ( x != "" ) x <- concat( header , x , footer )
        return(x)
    }
    prepend_indent <- function( x , n=1 ) {
        # adds n indents to start of each line after first
        return( gsub("\n",concat("\n",inden(n)),x,fixed=TRUE) )
    }

    get_symbol_dims <- function( symbol_name ) {
        return( symbols[[ symbol_name ]]$dims )
    }

    ###############################
    # reduce_sum hook BEGIN
    if ( threads > 1 ) {
        # first build list of data and parameters to pass to reduce_sum
        argument_names <- c()
        argument_types <- c()
        for ( i in length(symbols):1 ) {
            left_symbol <- names(symbols)[i]
            if ( symbols[[i]]$type %in% c("data","par","local") ) {
                if ( names(symbols)[i] != reduce_outcome ) {

                    # check for transpar
                    if ( symbols[[i]]$type=="local" ) {
                        if ( symbols[[i]]$block!="transpars" )
                            next
                    }

                    argument_names <- c( argument_names , names(symbols)[i] )
                    stan_type <- symbols[[i]]$dims

                    if ( class(stan_type)[1]=="character" )
                        stan_type <- list( stan_type , 1 )

                    if ( class(stan_type)[1]=="name" )
                        # symbol as dim, so likely a vector of varying effects
                        stan_type <- list( "vector" , 1 )

                    if ( stan_type[[1]]=="ordered" )
                        # ordered vector - pass to reduce_sum as vector
                        stan_type <- list( "vector" , 1 )

                    if ( stan_type[[1]]=="simplex" )
                        # simplex - pass to reduce_sum as vector
                        stan_type <- list( "vector" , 1 )

                    if ( stan_type[[1]]=="matrix" )
                        # matrix does not get [] on end
                        stan_type <- list( "matrix" , 1 )

                    if ( stan_type[[1]]=="cholesky_factor_corr" )
                        # simplex - pass to reduce_sum as vector
                        stan_type <- list( "matrix" , 1 )

                    if ( stan_type[[1]]!="vector" & stan_type[[2]] > 1 )
                         stan_type <- concat( stan_type[[1]] , "[]" )
                    else
                        stan_type <- stan_type[[1]]
                    argument_types <- c( argument_types , stan_type )
                }
            }
        }#i

        # now build code
        rsc <- concat( indent , "real reducer( \n" )
        reduce_outcome_type <- get_symbol_dims( reduce_outcome )[[1]]
        nn <- get_symbol_dims( reduce_outcome )[[2]]
        if ( nn > 1 & reduce_outcome_type != "vector" )
            reduce_outcome_type <- concat( reduce_outcome_type , "[]" )
        rsc <- concat( rsc , inden(3) , reduce_outcome_type , " " , reduce_outcome , ",\n" )
        rsc <- concat( rsc , inden(3) , "int start , int end , \n" )

        for ( i in 1:length(argument_names) ) {
            rsc <- concat( rsc , inden(3) , argument_types[i] , " " , argument_names[i] )
            if ( i < length(argument_names) ) rsc <- concat( rsc , ",\n" )
        }

        rsc <- concat( rsc , " ) { \n" )
        rsc <- concat( rsc , inden(1) , prepend_indent(m_f_reduce_declare,1) )

        # replace any XREDUCEOUTX replaceholder with outcome variable
        m_f_reduce <- gsub( "XREDUCEOUTX" , reduce_outcome , m_f_reduce , fixed=TRUE )

        rsc <- concat( rsc , prepend_indent(m_f_reduce,1) )
        rsc <- concat( rsc , "} \n" )
        m_funcs <- concat( m_funcs , rsc )

        # now add call in model block
        rsmc <- concat( inden(1) , "target += reduce_sum( reducer , " , reduce_outcome , " , " , floor(grain) , " , \n" )

        for ( i in 1:length(argument_names) ) {
            rsmc <- concat( rsmc , inden(3) , argument_names[i] )
            if ( i < length(argument_names) ) rsmc <- concat( rsmc , ",\n" )
        }

        rsmc <- concat( rsmc , " );\n" )
        model_body_append( rsmc , TRUE )
    }
    # reduce_sum hook END
    ###############################

    ###############################
    # compose it all
    m_funcs <- blockify( m_funcs , "functions{\n" , "}\n" )
    m_data <- blockify( m_data , "data{\n" , "}\n" )
    m_tdata1 <- blockify( m_tdata1 , "transformed data{\n" , "" )
    m_tdata2 <- blockify( m_tdata2 , "" , "}\n" )
    m_pars <- blockify( m_pars , "parameters{\n" , "}\n" )
    m_tpars1 <- blockify( m_tpars1 , "transformed parameters{\n" , "" )
    m_tpars2 <- blockify( m_tpars2 , "" , "}\n" )
    m_gq1 <- blockify( m_gq1 , "generated quantities{\n" , "" )
    m_gq2 <- blockify( m_gq2 , "" , "}\n" )
    
    model_code <- concat( 
        m_funcs ,
        m_data , 
        m_tdata1 , 
        m_tdata2 , 
        m_pars , 
        m_tpars1 , 
        m_tpars2 , 
        "model{\n" ,  
        m_model_declare , 
        m_model_txt , 
        "}\n" , 
        m_gq1 , 
        m_gq2 ,
        "\n" )

    stanfit <- NULL
    # use_pars <- c() --- declared earlier
    if ( missing(pars) ) {
        # return by default all pars, but no locals that might be in gen quants
        for ( i in 1:length(symbols) ) {
            if ( symbols[[i]]$type=="par" ) use_pars <- c( use_pars , names(symbols)[i] )
        }
    } else {
        use_pars <- pars
    }
    if ( log_lik==TRUE ) use_pars <- c( use_pars , "log_lik" )
    if ( length(pars_omit)>0 ) {
        # remove these names from use_pars
        idx <- which( pars_omit %in% use_pars )
        if ( length(idx)>0 ) use_pars <- use_pars[ -idx ]
    }

    # reverse order of use_pars, so names in formula order
    use_pars <- use_pars[ length(use_pars):1 ]

    cmdstanr_model_write <- function( the_model ) {
        # make temp name from model code md5 hash
        require( digest , quietly=TRUE )
        file_patt <- file.path( tempdir() , concat("ulam_cmdstanr_",digest(the_model,"md5")) )
        #file <- tempfile("ulam_cmdstanr",fileext=".stan")
        file_stan <- concat( file_patt , ".stan" )
        fileConn <- file( file_stan )
        writeLines( the_model , fileConn )
        close(fileConn)
        file_exe <- character()
        do_compile <- TRUE
        if ( file.exists(file_patt) ) {
            file_exe <- file_patt
            do_compile <- TRUE # force for now because of odd bug
        }
        return(list(file_stan,file_exe,do_compile))
    }

    #if ( threads>1 ) 
    cpp_options[['stan_threads']] <- TRUE

    # dangerous compile settings that can improve run times
    if ( cpp_fast==TRUE ) {
        cpp_options[['STAN_NO_RANGE_CHECKS']] <- TRUE
        cpp_options[['STAN_CPP_OPTIMS']] <- TRUE
    }

    # fire lasers pew pew
    if ( sample==TRUE ) {
        if ( length(start)==0 ) {
            # without explicit start list
            if ( prev_stanfit==FALSE ) {
                if ( cmdstan==FALSE )
                    # rstan interface
                    stanfit <- stan( model_code = model_code , data = data , pars=use_pars , chains=chains , cores=cores , iter=iter , control=control , warmup=warmup , ... )
                else {
                    # use cmdstanr interface
                    require( cmdstanr , quietly=TRUE )
                    filex <- cmdstanr_model_write( model_code )
                    mod <- cmdstan_model(
                        stan_file=filex[[1]],
                      # exe_file=filex[[2]],
                        compile=filex[[3]],
                        cpp_options=cpp_options )
                    # set_num_threads( threads )
                    # iter means only post-warmup samples for cmdstanr
                    # so need to compute iter explicitly
                    cmdstanfit <- mod$sample( 
                        data=data , 
                        chains=chains , parallel_chains=cores , 
                        iter_warmup=warmup, iter_sampling=floor(iter-warmup) , 
                        save_warmup=TRUE ,
                        adapt_delta=as.numeric(control[['adapt_delta']]) ,
                        threads_per_chain=threads , ... )
                    if ( rstanout==TRUE )
                        stanfit <- rstan::read_stan_csv(cmdstanfit$output_files())
                    else
                        stanfit <- cmdstanfit
                }
            } else {
                if ( cmdstan==FALSE ) {
                    # rstan interface, previous stanfit object
                    stanfit <- stan( fit = prev_stanfit_object , data = data , pars=use_pars , 
                         chains=chains , cores=cores , iter=iter , control=control , warmup=warmup , ... )
                } else {
                    # SAME AS ABOVE FOR NOW - how to referece exe?
                    # use cmdstanr interface
                    require( cmdstanr , quietly=TRUE )
                    filex <- cmdstanr_model_write( model_code )
                    mod <- cmdstan_model(
                        stan_file=filex[[1]],
                      # exe_file=filex[[2]],
                        compile=filex[[3]],
                        cpp_options=cpp_options )
                    # set_num_threads( threads )
                    # iter means only post-warmup samples for cmdstanr
                    # so need to compute iter explicitly
                    cmdstanfit <- mod$sample( 
                        data=data , 
                        chains=chains , parallel_chains=cores , 
                        iter_warmup=warmup, iter_sampling=floor(iter-warmup) , 
                        save_warmup=TRUE ,
                        adapt_delta=as.numeric(control[['adapt_delta']]) ,
                        threads_per_chain=threads , ... )
                    if ( rstanout==TRUE )
                        stanfit <- rstan::read_stan_csv(cmdstanfit$output_files())
                    else
                        stanfit <- cmdstanfit
                }
            }
        } else {
            # WITH explicit start list
            # Still needs cmdstanr interface
            f_init <- "random"
            if ( class(start)=="list" ) f_init <- function() return(start)
            if ( class(start)=="function" ) f_init <- start
            if ( prev_stanfit==FALSE )
                stanfit <- stan( model_code = model_code , data = data , pars=use_pars , 
                         chains=chains , cores=cores , iter=iter , control=control , init=f_init , warmup=warmup , ... )
            else
                stanfit <- stan( fit = prev_stanfit_object , data = data , pars=use_pars , 
                         chains=chains , cores=cores , iter=iter , control=control , init=f_init , warmup=warmup , ... )
        }
    }

    formula_parsed <- list(
                formula = flist,
                symbols = symbols,
                indices = indices,
                link_funcs = link_funcs,
                symbol_graph = symbol_graph,
                impute_bank = impute_bank,
                marginalize_bank = marginalize_bank
            )

    if ( sample==FALSE ) {
        result <- list(
            #stanfit = stanfit,
            formula = flist.orig,
            model = model_code,
            data = data,
            formula_parsed = formula_parsed
        )
    } else {

        # compute coef vector and vcov
        if ( rstanout==TRUE | cmdstan==FALSE ) {
            s <- summary(stanfit)$summary
            s <- s[ -which( rownames(s)=="lp__" ) , ] # clean out lp__
            if ( log_lik==TRUE ) {
                # clean out the log_lik entries - possibly very many
                idx <- grep( "log_lik[" , rownames(s) , fixed=TRUE )
                s <- s[ -idx , ]
            }
            if ( !is.null(dim(s)) ) {
                coef <- s[,1]
                varcov <- matrix( s[,3]^2 , ncol=1 )
                rownames(varcov) <- rownames(s)
            } else {
                coef <- s[1]
                varcov <- matrix( s[3]^2 , ncol=1 )
                #names(coef) <- names(start[[1]])
            }
        } else {
            coef <- 1:2
            varcov <- matrix(NA,2,2)
        }

        result <- new( "ulam" , 
                call = match.call(), 
                model = model_code,
                #stanfit = stanfit,
                coef = coef,
                vcov = varcov,
                data = data,
                start = list( start ),
                pars = use_pars,
                formula = flist.orig,
                formula_parsed = formula_parsed
            )
        if ( rstanout==TRUE | cmdstan==FALSE ) {
            attr(result,"stanfit") <- stanfit
        } else {
            attr(result,"cstanfit") <- stanfit
        }
        
        attr(result,"generation") <- "ulam2018"
        if ( nobs_save > 0 ) attr(result,"nobs") <- nobs_save

    }

    # check for file argument
    if ( !is.null(file) ) {
        rds_file_name <- concat( file , ".rds" )
        if ( !file.exists(rds_file_name) ) {
            # save it
            if ( messages==TRUE ) message(concat("Saving result as ",rds_file_name))
            saveRDS( result , file=rds_file_name )
        }
    }

    return( result )

}

# new link function
setMethod( "link" , "ulam" , function(fit,...) link_ulam(fit,...) )
link_ulam <- function( fit , data , post , simplify=TRUE , symbols , ... ) {

    nest_commas <- function( f , symbols ) {
        if ( class( f )=="call" || class( f )=="(" ) {
            if ( as.character( f[[1]] )=="[" ) {
                # bracket call, so maybe insert comma (TRUE)
                if ( as.character( f[[2]] ) %in% symbols ) {
                    the_call <- as.list( f )
                    f <- as.call( unlist( list( the_call[1:2] , quote(TRUE) , the_call[3:length(the_call)] ) , recursive=FALSE ) )
                }
            } else {
                # need to drill down
                for ( i in 1:length(f) )
                    f[[i]] <- nest_commas( f[[i]] , symbols )
            }
        }
        return( f )
    }

    use_orig_data <- FALSE
    if ( missing(data) ) {
        data <- fit@data
        use_orig_data <- TRUE
    }
    if ( missing(post) ) post <- extract.samples(fit@stanfit)

    # how many link functions?
    if ( !is.null(fit@formula_parsed$link_funcs) ) {
        n_links <- length( fit@formula_parsed$link_funcs )
    } else {
        stop( "No link functions found. Maybe there were no linear models in the formula list?" )
    }

    out <- list()

    if ( missing(symbols) ) {
        # use them all
        symbols <- names( fit@formula_parsed$link_funcs )
    }

    for ( j in 1:n_links ) {

        if ( !(names( fit@formula_parsed$link_funcs )[j] %in% symbols) ) next

        # get number of cases for this symbol
        if ( use_orig_data==TRUE )
            n_cases <- fit@formula_parsed$link_funcs[[ j ]]$N
        else
            # guess from length of data
            n_cases <- length(data[[1]])

        # check whether this function contains symbols that need some dimensions added for samples
        symbols_so_far <- names(out)
        f <- fit@formula_parsed$link_funcs[[ j ]]$func
        for ( jj in symbols_so_far ) {
            # more than one dim? if so, then may need extra dim when referenced
            if ( length(dim( out[[jj]] ))>1 ) {
                f <- nest_commas( f , jj )
            }
        }

        l <- sapply( 1:n_cases , function(i) {
                the_env <- post
                #the_env <- unlist( list( the_env , data[i,,drop=FALSE] ) , recursive=FALSE )
                the_env <- unlist( list( the_env , data , i=i ) , recursive=FALSE )
                if ( j > 1 ) {
                    # insert previous linear model results, in case used in later ones
                    for ( k in 1:length(out) ) {
                        lv_name <- names(out)[k]
                        #the_env[[ lv_name ]] <- out[[k]][ , i ] # all samples, i-th case
                        the_env[[ lv_name ]] <- out[[k]]
                    }
                }
                eval( f , envir=the_env ) 
            })

        # an inverse link to apply?
        if ( !is.null( fit@formula_parsed$link_funcs[[ j ]]$inv_link ) ) {
            l <- do.call( fit@formula_parsed$link_funcs[[ j ]]$inv_link , list(l) )
        }

        out[[ names( fit@formula_parsed$link_funcs )[j] ]] <- l

    }#j

    n_links <- length(symbols)

    if ( simplify==TRUE && n_links==1 ) out <- out[[1]]

    return(out)
}

# sim method assumes outcome var is first var on left in formula line 1
# can use variable argument to specify any other variable
setMethod( "sim" , "ulam" , function(fit,...) sim_ulam(fit,...) )
sim_ulam <- function( fit , data , post , variable , n=1000 , replace=list() , ... ) {
    
    ########################################
    # check arguments
    if ( missing(data) ) {
        data <- fit@data
    }

    if ( n==0 ) {
        n <- stan_total_samples(fit@stanfit)
    } else {
        tot_samples <- stan_total_samples(fit@stanfit)
        n <- min(n,tot_samples)
    }
    
    if ( missing(post) ) 
        post <- extract.samples(fit,n=n)

    # get linear model values from link
    # use our posterior samples, so later parameters have right correlation structure with link values
    # don't flatten result, so we end up with a named list, even if only one element
    pred <- link_ulam( fit , data=data , post=post , simplify=FALSE , ... )
    
    # extract likelihood, assuming it is first element of formula
    if ( missing(variable) ) {
        lik <- fit@formula[[1]]
        # discover outcome
        outcome <- as.character(lik[[2]])
        # discover likelihood function
        flik <- as.character(lik[[3]][[1]])
    } else {
        # named outcome, so find it
        for ( i in 1:length(fit@formula) ) {
            outcome <- as.character( fit@formula[[i]][[2]] )
            if ( outcome==variable ) {
                lik <- fit@formula[[i]]
                flik <- as.character( lik[[3]][[1]] )
                break
            }
        }#i
    }

    # check whether we must convert from Stan-name distribution to R-name distribution
    first_char <- substr( flik , 1 , 1 )
    if ( first_char != "d" ) {
        # loop over templates and get matching dfoo name
        for ( ii in 1:length( ulam_dists ) ) {
            aStanName <- ulam_dists[[ii]]$Stan_name
            if ( aStanName==flik ) {
                flik <- ulam_dists[[ii]]$R_name
                break
            }
        }#ii
    }

    # get simulation partner function
    rlik <- flik
    substr( rlik , 1 , 1 ) <- "r"

    # pull out parameters in likelihood
    pars <- vector(mode="list",length=length(lik[[3]])-1)
    for ( i in 1:length(pars) ) {
        pars[[i]] <- lik[[3]][[i+1]]
    }
    pars <- paste( pars , collapse=" , " )
    # build expression to evaluate
    n_cases <- length(data[[1]])
    xeval <- paste( rlik , "(" , n_cases , "," , pars , ")" , collapse="" )
    
    # simulate outcomes
    sim_out <- matrix( NA , nrow=n , ncol=n_cases )
    # need to select out s-th sample for each parameter in post
    # only need top-level parameters, so can coerce to data.frame?
    post2 <- as.data.frame(post)
    for ( s in 1:n ) {
        # build environment
        #ndims <- length(dim(pred[[1]]))
        #if ( ndims==2 )
        #    elm <- pred[[1]][s,] # extract s-th sample case calculations
        #if ( ndims==3 )
        #    elm <- pred[[1]][s,,]
        #elm <- list(elm)
        #names(elm)[1] <- names(pred)[1]
        
        elm <- list()
        for ( j in 1:length(pred) ) {
            ndims <- length(dim(pred[[j]]))
            if ( ndims==2 )
                elm[[j]] <- pred[[j]][s,]
            if ( ndims==3 )
                elm[[j]] <- pred[[j]][s,,]
        }
        names(elm) <- names(pred)
        
        e <- list( as.list(data) , as.list(post2[s,]) , as.list(elm) )
        e <- unlist( e , recursive=FALSE )
        # evaluate
        if ( flik=="dordlogit" ) {
            n_outcome_vals <- dim( pred[[1]] )[3]
            probs <- pred[[1]][s,,]
            sim_out[s,] <- sapply( 
                1:n_cases ,
                function(i)
                    sample( 1:n_outcome_vals , size=1 , replace=TRUE , prob=probs[i,] )
            )
        } else {
            sim_out[s,] <- eval(parse(text=xeval),envir=e)
        }
    }
    
    return(sim_out)
}



