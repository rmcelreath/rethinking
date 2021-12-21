# functions for working with cmdstan and cmdstanr

# write model code as character vector to file, so cmdstan_model can read it
cmdstanr_model_write <- function( the_model ) {
        # make temp name from model code md5 hash
        require( digest , quietly=TRUE )
        file_patt <- file.path( tempdir() , concat("rt_cmdstanr_",digest(the_model,"md5")) )
        #file <- tempfile("ulam_cmdstanr",fileext=".stan")
        file_stan <- concat( file_patt , ".stan" )
        fileConn <- file( file_stan )
        writeLines( the_model , fileConn )
        close(fileConn)
        return(file_stan)
    }

# wrapper to fit stan model using cmdstanr and return rstan stanfit object
cstan <- function( file , model_code , data=list() , chains=1 , cores=1 , iter=1000 , warmup , threads=1 , control=list(adapt_delta=0.95) , cpp_options=list() , save_warmup=TRUE , cpp_fast=FALSE , rstan_out=TRUE , pars , ... ) {

    if ( threads>1 ) cpp_options[['stan_threads']] <- TRUE

    # dangerous compile settings that can improve run times
    if ( cpp_fast==TRUE ) {
        cpp_options[['STAN_NO_RANGE_CHECKS']] <- TRUE
        cpp_options[['STAN_CPP_OPTIMS']] <- TRUE
    }

    if ( missing(file) & !missing(model_code) ) {
        file <- cmdstanr_model_write( model_code )
    }

    require(cmdstanr,quietly=TRUE)
    mod <- cmdstan_model( file , compile=TRUE , cpp_options=cpp_options )

    if ( missing(warmup) ) {
        samp <- floor(iter/2)
        warm <- floor(iter/2)
    } else {
        samp <- iter - warmup
        warm <- warmup
    } 

    # pull out any control arguments
    carg_adapt_delta <- 0.95
    if ( !is.null( control[['adapt_delta']] ) )
        carg_adapt_delta <- as.numeric(control[['adapt_delta']])
    carg_max_treedepth <- 11
    if ( !is.null( control[['max_treedepth']] ) )
        carg_max_treedepth <- as.numeric(control[['max_treedepth']])

    # sample
    if ( threads > 1 )
        cmdstanfit <- mod$sample( data=data , 
            chains=chains , 
            parallel_chains=cores , 
            iter_sampling=samp , iter_warmup=warm , 
            adapt_delta=carg_adapt_delta , 
            max_treedepth=carg_max_treedepth , 
            threads_per_chain=threads ,
            save_warmup=save_warmup , ... )
    else
        cmdstanfit <- mod$sample( data=data , 
            chains=chains , 
            parallel_chains=cores , 
            iter_sampling=samp , iter_warmup=warm , 
            adapt_delta=carg_adapt_delta , 
            max_treedepth=carg_max_treedepth ,
            save_warmup=save_warmup , ... )

    # coerce to stanfit object
    if ( rstan_out==TRUE )
        return( rstan::read_stan_csv(cmdstanfit$output_files()) )
    else
        return(cmdstanfit)
}

# override rstan's stan() function to use cmdstanr?

stan <- function( ... ) {
    if ( ulam_options$use_cmdstan==TRUE ) {
        cstan( ... )
    } else {
        rstan::stan( ... )
    }
}

# in place tests
if ( FALSE ) {

    library(rethinking)

    the_code <- "
    data{
        int N;
        real y[N];
    }
    parameters{
        real mu;
    }
    model{
        mu ~ normal(0,1);
        y ~ normal(mu,1);
    }
    "

    m <- cstan( model_code=the_code , data=list(N=2,y=rnorm(2)) , chains=4 , rstan_out=FALSE )

    m$summary( variables="mu" , "mean" , "sd" , ~quantile(.x, probs = c(0.4, 0.6)) , "rhat" , "ess_bulk" )

    precis(m)

}
