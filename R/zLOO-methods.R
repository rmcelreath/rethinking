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
        looIC <- as.vector( loo_list$pointwise[,4] )
        lppd <- as.vector( loo_list$pointwise[,1] )
        pD <- as.vector( loo_list$pointwise[,3] )
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
        ll_matrix <- loo::extract_log_lik( object , merge_chains=FALSE )
    else
        stop(concat("Log-likelihood matrix '",log_lik,"'' not found."))

    rel_n_eff <- loo::relative_eff(exp(ll_matrix))

    loo_list <- loo::loo( ll_matrix , r_eff = rel_n_eff )

    if ( pointwise==TRUE ) {
        looIC <- as.vector( loo_list$pointwise[,4] )
        lppd <- as.vector( loo_list$pointwise[,1] )
        pD <- as.vector( loo_list$pointwise[,3] )
    } else {
        looIC <- loo_list$estimates[3,1]
        lppd <- loo_list$estimates[1,1]
        pD <- loo_list$estimates[2,1]
    }
    attr(looIC,"lppd") = lppd
    attr(looIC,"pLOO") = pD
    
    n_tot <- ncol(ll_matrix)
    #attr(looIC,"se") = try(sqrt( n_tot*var(as.vector( loo_list$pointwise[,4] )) ))
    attr(looIC,"se") = loo_list$estimates[3,2]
    
    return(looIC)
})

setMethod("LOO", "ulam",
function( object , n=0 , refresh=0.1 , pointwise=FALSE , log_lik="log_lik" , ... ) {
    return( LOO(object@stanfit,n=n,refresh=refresh,pointwise=pointwise,log_lik=log_lik,...) )
})

# extracts log_lik matrix from samples in a list
setMethod("LOO", "list",
function( object , n=0 , refresh=0.1 , pointwise=FALSE , log_lik="log_lik" , ... ) {
    if ( !is.null(object[[log_lik]]) )
        ll_matrix <- object[[log_lik]]
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
function( object , n=1000 , refresh=0 , pointwise=FALSE , ... ) {
    
    # get log-likelihood matrix from sim method
    loglik_matrix <- sim(object,n=n,refresh=refresh,ll=TRUE,...)
    
    loo_list <- suppressWarnings( loo::loo(loglik_matrix) )
    
    if ( pointwise==TRUE ) {
        looIC <- as.vector( loo_list$pointwise[,4] )
        lppd <- as.vector( loo_list$pointwise[,1] )
        pD <- as.vector( loo_list$pointwise[,3] )
    } else {
        looIC <- loo_list$estimates[3,1]
        lppd <- loo_list$estimates[1,1]
        pD <- loo_list$estimates[2,1]
    }
    attr(looIC,"lppd") = lppd
    attr(looIC,"pLOO") = pD
    attr(looIC,"diagnostics") = loo_list$diagnostics
    
    n_tot <- ncol(loglik_matrix)
    #attr(looIC,"se") = try(sqrt( n_tot*var2(as.vector( loo_list$pointwise[,4] )) ))
    attr(looIC,"se") = loo_list$estimates[3,2]
    
    return(looIC)

})

