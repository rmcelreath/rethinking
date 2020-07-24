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
cstan <- function( file , model_code , data=list() , chains=1 , cores=1 , iter=1000 , threads=1 , control=list(adapt_delta=0.95) , cpp_options=list() , ... ) {

    if ( threads>1 ) cpp_options[['stan_threads']] <- TRUE

    if ( missing(file) & !missing(model_code) ) {
        file <- cmdstanr_model_write( model_code )
    }

    require(cmdstanr,quietly=TRUE)
    mod <- cmdstan_model( file , compile=TRUE , cpp_options=cpp_options )

    if ( threads > 1 )
        cmdstanfit <- mod$sample( data=data , chains=chains , parallel_chains=cores , iter_sampling=floor(iter/2) , adapt_delta=as.numeric(control[['adapt_delta']]) , threads_per_chain=threads , ... )
    else
        cmdstanfit <- mod$sample( data=data , chains=chains , parallel_chains=cores , iter_sampling=floor(iter/2) , adapt_delta=as.numeric(control[['adapt_delta']]) , ... )

    stanfit <- rstan::read_stan_csv(cmdstanfit$output_files())

    return(stanfit)
}

