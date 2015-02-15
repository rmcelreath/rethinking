# build ensemble of samples using DIC/WAIC weights
ensemble <- function( ... , data , n=1e3 , WAIC=TRUE , refresh=0 ) {
    # retrieve list of models
    L <- list(...)
    if ( is.list(L[[1]]) && length(L)==1 )
        L <- L[[1]]
    # retrieve model names from function call
    mnames <- match.call()
    mnames <- as.character(mnames)[2:(length(L)+1)]
    if ( length(L)>1 ) {
        ictab <- compare( ... , WAIC=WAIC , refresh=refresh , n=n , sort=FALSE )
        rownames(ictab@output) <- mnames
        weights <- ictab@output$weight
    } else {
        ictab <- NA
        weights <- 1
    }
    
    # compute number of predictions per model
    idx <- round( weights * n )
    
    # generate simulated predictions for each model
    # then compose ensemble using weights
    if ( missing(data) ) {
        link.list <- lapply( L , function(m) link(m , n=n , refresh=refresh ) )
        sim.list <- lapply( L , function(m) sim(m , n=n , refresh=refresh ) )
    } else {
        link.list <- lapply( L , function(m) link(m , data=data , n=n , refresh=refresh ) )
        sim.list <- lapply( L , function(m) sim(m , data=data , n=n , refresh=refresh ) )
    }
    #print(str(sim.list))
    
    # compose weighted predictions
    # calculate indices corresponding to each model
    idx_start <- rep(1,length(idx))
    idx_end <- idx
    if ( length(L)>1 )
        for ( i in 2:length(idx) ) {
            idx_start[i] <- min( idx_start[i-1] + idx[i-1] , n )
            idx_end[i] <- min( idx_end[i-1] + idx[i] , n )
        }
    #print(idx)
    #print(idx_start)
    #print(idx_end)
    link_out <- link.list[[1]]
    sim_out <- sim.list[[1]]
    for ( i in 1:length(idx) ) {
        if ( idx[i]>0 ) {
            idxrange <- idx_start[i]:idx_end[i]
            link_out[idxrange,] <- link.list[[i]][idxrange,]
            sim_out[idxrange,] <- sim.list[[i]][idxrange,]
        }
    }
    result <- list( link=link_out , sim=sim_out )
    names(weights) <- mnames
    idxtab <- cbind( idx_start , idx_end )
    rownames(idxtab) <- mnames
    attr(result,"weights") <- weights
    attr(result,"indices") <- idxtab
    attr(result,"ictab") <- ictab
    result
}

