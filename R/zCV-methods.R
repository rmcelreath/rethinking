## methods for explicit cross-validation

dims <- function(x) {
    z <- dim(x)
    if ( is.null(z) ) return(1) # must be 1D vector
    return( length(z) )
}

cv_data_slice <- function( do , N , exclude , include ) {
    dlo <- list()
    if ( missing(include) ) include <- (1:N)[-exclude]
    for ( j in 1:length(do) ) {
        vname <- names(do)[j]
        if ( dims(do[[j]])==1 ) {
            if ( length(do[[j]])==N ) 
                dlo[[vname]] <- do[[vname]][ include ]
            else
                dlo[[vname]] <- do[[vname]] # all of it
        } else {
            # matrix?
            # check if rows of length N
            the_dims <- dim( do[[j]] )
            if ( the_dims[1]==N ) {
                dlo[[vname]] <- do[[vname]][ include , ]
                # catch collapsing matrix problem
                if ( class(dlo[[vname]])=="numeric" )
                    dlo[[vname]] <- matrix( dlo[[vname]] , nrow=length(include) , ncol=ncol(do[[vname]]) )
            } else {
                # not sure what to do
                dlo[[vname]] <- do[[vname]]
            }
        }
    }#j
    return(dlo)
}

cv_quap <- function( quap_model , lno=1 , pw=FALSE , cores=1 , ... ) {
    outcome <- deparse( quap_model@formula[[1]][[2]] )
    N <- length( quap_model@data[[outcome]] )
    # number of models to fit
    m <- floor( N / lno )
    mlist <- list()
    do <- quap_model@data # could be data.frame or list
    d_list <- list()
    for ( i in 1:m ) {
        leave_out <- ( 1 + (i-1)*lno ):( lno + (i-1)*lno )
        d_list[[i]] <- cv_data_slice( do , N , exclude=leave_out )
    }
    mlist <- mclapply( 1:m , function(i) quap( quap_model@formula , data=d_list[[i]] , ... ) , mc.cores=cores , mc.allow.recursive=TRUE )
    # now for each model, predict the left out shard
    lppd_list <- list()
    leave_out <- list()
    for ( i in 1:m ) {
        ix <- ( 1 + (i-1)*lno ):( lno + (i-1)*lno )
        leave_out[[i]] <- cv_data_slice( do , N , include=ix )
        #lppd_list[[i]] <- lppd( mlist[[i]] , data=cv_data_slice( do , N , include=leave_out ) )
    }
    lppd_list <- mclapply( 1:m , function(i) lppd( mlist[[i]] , data=leave_out[[i]] , mc.cores=cores ,  mc.allow.recursive=TRUE ) )
    if ( lno==1 ) lppd_list <- unlist( lppd_list )
    if ( pw==FALSE ) lppd_list <- sum( unlist(lppd_list) )
    # result
    return( lppd_list )
}

