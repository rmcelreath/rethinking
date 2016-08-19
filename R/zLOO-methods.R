###################################
# compute Vehtari's PSIS-LOO from map or map2stan fit

setGeneric("LOO",
function( object , n=1000 , refresh=0.1 , pointwise=FALSE , ... ) {
    message( concat("No LOO method for object of class '",class(object),"'. Returning AIC instead.") )
    AIC(object)
}
)

setMethod("LOO", "map2stan",
function( object , n=0 , refresh=0.1 , pointwise=FALSE , ... ) {
    
    # get log-likelihood matrix from WAIC method
    loglik_matrix <- WAIC(object,n=n,refresh=refresh,loglik=TRUE,...)
    
    loo_list <- loo::loo(loglik_matrix)

    if ( pointwise==TRUE ) {
        looIC <- as.vector( loo_list$pointwise[,3] )
        lppd <- as.vector( loo_list$pointwise[,1] )
        pD <- as.vector( loo_list$pointwise[,2] )
    } else {
        looIC <- loo_list$looic
        lppd <- loo_list$elpd_loo
        pD <- loo_list$p_loo
    }
    attr(looIC,"lppd") = lppd
    attr(looIC,"pLOO") = pD
    
    n_tot <- ncol(loglik_matrix)
    attr(looIC,"se") = try(sqrt( n_tot*var2(as.vector( loo_list$pointwise[,3] )) ))
    
    return(looIC)

})

# extracts log_lik matrix from stanfit and computes LOO
setMethod("LOO", "stanfit",
function( object , n=0 , refresh=0.1 , pointwise=FALSE , log_lik="log_lik" , ... ) {
    if ( !is.null(extract.samples(object)[[log_lik]]) )
        ll_matrix <- extract.samples(object)[[log_lik]]
    else
        stop(concat("Log-likelihood matrix '",log_lik,"'' not found."))

    loo_list <- loo::loo(ll_matrix)

    if ( pointwise==TRUE ) {
        looIC <- as.vector( loo_list$pointwise[,3] )
        lppd <- as.vector( loo_list$pointwise[,1] )
        pD <- as.vector( loo_list$pointwise[,2] )
    } else {
        looIC <- loo_list$looic
        lppd <- loo_list$elpd_loo
        pD <- loo_list$p_loo
    }
    attr(looIC,"lppd") = lppd
    attr(looIC,"pLOO") = pD
    
    n_tot <- ncol(ll_matrix)
    attr(looIC,"se") = try(sqrt( n_tot*var2(as.vector( loo_list$pointwise[,3] )) ))
    
    return(looIC)
})

setMethod("LOO", "map",
function( object , n=1000 , refresh=0.1 , pointwise=FALSE , ... ) {
    
    # get log-likelihood matrix from sim method
    loglik_matrix <- sim(object,n=n,refresh=refresh,ll=TRUE,...)
    
    loo_list <- loo::loo(loglik_matrix)
    
    if ( pointwise==TRUE ) {
        looIC <- as.vector( loo_list$pointwise[,3] )
        lppd <- as.vector( loo_list$pointwise[,1] )
        pD <- as.vector( loo_list$pointwise[,2] )
    } else {
        looIC <- loo_list$looic
        lppd <- loo_list$elpd_loo
        pD <- loo_list$p_loo
    }
    attr(looIC,"lppd") = lppd
    attr(looIC,"pLOO") = pD
    
    n_tot <- ncol(loglik_matrix)
    attr(looIC,"se") = try(sqrt( n_tot*var2(as.vector( loo_list$pointwise[,3] )) ))
    
    return(looIC)

})

