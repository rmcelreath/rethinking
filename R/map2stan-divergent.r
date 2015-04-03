# extracts n_divergent from stan fit
divergent <- function( fit , warmup=FALSE ) {
    if ( class(fit)=="map2stan" ) fit <- fit@stanfit
    x <- rstan::get_sampler_params(fit)
    if ( warmup==FALSE ) {
        nwarmup <- fit@stan_args[[1]]$warmup
        niter <- fit@stan_args[[1]]$iter
        n <- sapply( x , function(ch) sum(ch[(nwarmup+1):niter,5]) )
    } else {
        n <- sapply( x , function(ch) sum(ch[,5]) )
    }
    sum(n)
}
