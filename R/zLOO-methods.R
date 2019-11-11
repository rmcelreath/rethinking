###################################
# compute Vehtari's PSIS-LOO from map or map2stan fit

setGeneric("LOO",
function( object , n=1000 , refresh=0.1 , pointwise=FALSE , ... ) {
    message( concat("No PSIS method for object of class '",class(object),"'. Returning AIC instead.") )
    AIC(object)
}
)

setGeneric("PSIS",
function( object , n=1000 , refresh=0.1 , pointwise=FALSE , ... ) {
    message( concat("No PSIS method for object of class '",class(object),"'. Returning AIC instead.") )
    AIC(object)
}
)

PSIS_map2stan <- function( object , n=0 , refresh=-1 , pointwise=FALSE , warn=TRUE , ... ) {
    
    # get log-likelihood matrix from WAIC method
    ll_matrix <- WAIC(object,n=n,refresh=refresh,loglik=TRUE,...)
    
    rel_n_eff <- loo::relative_eff( exp(ll_matrix) , chain_id=rep(1,nrow(ll_matrix)) )
    loo_list <- suppressWarnings( loo::loo( ll_matrix , r_eff = rel_n_eff ) )

    #loo_list <- loo::loo(loglik_matrix)

    if ( pointwise==TRUE ) {
        looIC <- as.vector( loo_list$pointwise[,4] )
        lppd <- as.vector( loo_list$pointwise[,1] )
        pD <- as.vector( loo_list$pointwise[,3] )
    } else {
        looIC <- loo_list$estimates[3,1]
        lppd <- loo_list$estimates[1,1]
        pD <- loo_list$estimates[2,1]
    }
    #attr(looIC,"lppd") = lppd
    #attr(looIC,"pLOO") = pD
    #attr(looIC,"diagnostics") = loo_list$diagnostics

    if ( warn==TRUE ) rethinking:::xcheckLOOk( loo_list$diagnostics$pareto_k )
    
    n_tot <- ncol(ll_matrix)
    #attr(looIC,"se") = try(sqrt( n_tot*var2(as.vector( loo_list$pointwise[,3] )) ))

    result <- data.frame( 
        PSIS=looIC , 
        lppd=lppd , 
        penalty=pD ,
        std_err=loo_list$estimates[3,2] )

    if ( pointwise==TRUE ) {
        result$k = loo_list$diagnostics$pareto_k
    }
    
    return(result)

}
setMethod("LOO","map2stan", PSIS_map2stan )
setMethod("PSIS", "map2stan", PSIS_map2stan )

# extracts log_lik matrix from stanfit and computes LOO
PSIS_stanfit <- function( object , n=0 , refresh=0.1 , pointwise=FALSE , log_lik="log_lik" , warn=TRUE , ... ) {
    if ( !is.null(extract.samples(object)[[log_lik]]) )
        ll_matrix <- loo::extract_log_lik( object , merge_chains=FALSE )
    else
        stop(concat("Log-likelihood matrix '",log_lik,"'' not found."))

    rel_n_eff <- loo::relative_eff(exp(ll_matrix))

    loo_list <- suppressWarnings( loo::loo( ll_matrix , r_eff = rel_n_eff ) )

    if ( pointwise==TRUE ) {
        looIC <- as.vector( loo_list$pointwise[,4] )
        lppd <- as.vector( loo_list$pointwise[,1] )
        pD <- as.vector( loo_list$pointwise[,3] )
    } else {
        looIC <- loo_list$estimates[3,1]
        lppd <- loo_list$estimates[1,1]
        pD <- loo_list$estimates[2,1]
    }
    #attr(looIC,"lppd") = lppd
    #attr(looIC,"pLOO") = pD
    #attr(looIC,"diagnostics") = loo_list$diagnostics

    if ( warn==TRUE ) rethinking:::xcheckLOOk( loo_list$diagnostics$pareto_k )
    
    n_tot <- ncol(ll_matrix)
    #attr(looIC,"se") = try(sqrt( n_tot*var(as.vector( loo_list$pointwise[,4] )) ))
    #attr(looIC,"se") = loo_list$estimates[3,2]

    result <- data.frame( 
        PSIS=looIC , 
        lppd=lppd , 
        penalty=pD ,
        std_err=loo_list$estimates[3,2] )

    if ( pointwise==TRUE ) {
        result$k = loo_list$diagnostics$pareto_k
    }
    
    return(result)
}
setMethod("PSIS", "stanfit", PSIS_stanfit )
setMethod("LOO", "stanfit", PSIS_stanfit )

PSIS_ulam <- function( object , n=0 , refresh=0.1 , pointwise=FALSE , log_lik="log_lik" , warn=TRUE , ... ) {
    return( PSIS(object@stanfit,n=n,refresh=refresh,pointwise=pointwise,log_lik=log_lik,warn=warn,...) )
}
setMethod("LOO", "ulam", PSIS_ulam)
setMethod("PSIS", "ulam", PSIS_ulam)

# extracts log_lik matrix from samples in a list
PSIS_list <- function( object , n=0 , refresh=0.1 , pointwise=FALSE , log_lik="log_lik" , warn=TRUE , ... ) {
    if ( !is.null(object[[log_lik]]) )
        ll_matrix <- object[[log_lik]]
    else
        stop(concat("Log-likelihood matrix '",log_lik,"'' not found."))

    loo_list <- suppressWarnings( loo::loo(ll_matrix) )

    if ( pointwise==TRUE ) {
        looIC <- as.vector( loo_list$pointwise[,4] )
        lppd <- as.vector( loo_list$pointwise[,1] )
        pD <- as.vector( loo_list$pointwise[,3] )
    } else {
        looIC <- loo_list$estimates[3,1]
        lppd <- loo_list$estimates[1,1]
        pD <- loo_list$estimates[2,1]
    }
    #attr(looIC,"lppd") = lppd
    #attr(looIC,"pLOO") = pD
    #attr(looIC,"diagnostics") = loo_list$diagnostics

    if ( warn==TRUE ) rethinking:::xcheckLOOk( loo_list$diagnostics$pareto_k )
    
    n_tot <- ncol(ll_matrix)
    #attr(looIC,"se") = try(sqrt( n_tot*var2(as.vector( loo_list$pointwise[,3] )) ))

    result <- data.frame( 
        PSIS=looIC , 
        lppd=lppd , 
        penalty=pD ,
        std_err=loo_list$estimates[3,2] )

    if ( pointwise==TRUE ) {
        result$k = loo_list$diagnostics$pareto_k
    }
    
    return(result)
}
setMethod("LOO", "list",PSIS_list)
setMethod("PSIS", "list",PSIS_list)

PSIS_quap <- function( object , n=1000 , refresh=0 , pointwise=FALSE , warn=TRUE , ... ) {
    
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

    #attr(looIC,"lppd") = lppd
    #attr(looIC,"pLOO") = pD
    #attr(looIC,"diagnostics") = loo_list$diagnostics

    #object_name <- deparse( match.call()[[2]] )
    if ( warn==TRUE ) xcheckLOOk( loo_list$diagnostics$pareto_k )
    
    n_tot <- ncol(loglik_matrix)
    #attr(looIC,"se") = try(sqrt( n_tot*var2(as.vector( loo_list$pointwise[,4] )) ))

    result <- data.frame( 
        PSIS=looIC , 
        lppd=lppd , 
        penalty=pD ,
        std_err=loo_list$estimates[3,2] )

    if ( pointwise==TRUE ) {
        result$k = loo_list$diagnostics$pareto_k
    }

    return(result)

}
setMethod("LOO", "map", PSIS_quap)
setMethod("PSIS", "map", PSIS_quap)

# function to return pointwise k diagnostics
PSISk <- function( model , ... ) {
    round( PSIS(model,pointwise=TRUE)$k , 2 )
}

xcheckLOOk <- function( k_list ) {
    if ( any( k_list > 0.5 ) ) {
        if ( any( k_list > 1 ) ) {
            message( concat("Some Pareto k values are very high (>1). Set pointwise=TRUE to inspect individual points.") )
        } else {
            message( concat("Some Pareto k values are high (>0.5). Set pointwise=TRUE to inspect individual points.") )
        }
    }
}
