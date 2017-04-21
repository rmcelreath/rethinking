# build ensemble of samples using DIC/WAIC weights


#' Simulate an ensemble of posterior predictions
#'
#' Uses \code{link} and \code{sim} for a list of \code{map} or \code{map2stan}
#' model fits to construct Akaike weighted ensemble of predictions.
#'
#' This function calls \code{\link{link}} and \code{\link{sim}} for each fit
#' model given as input. The results are then combined into ensemble link and
#' simulation output, where samples from each model are represented in
#' proportion to the Akaike weights. Akaike weights are calculated by
#' \code{\link{compare}}, using \code{func}, unless an explicit vector
#' \code{weights} is provided. The values in \code{weights} will be normalized
#' to sum to one, if they do not already.
#'
#' Note that no averaging is done by this function. Parameters are not
#' averaged, and predictions are not averaged.
#'
#' @param ... \code{map} or \code{map2stan} models
#' @param data Optional data to compute predictions over, as in \code{link} and
#' \code{sim}
#' @param n Number of samples to draw from posterior for each model
#' @param func Function to use in computing criterion for model weights
#' @param weights Optional vector of weights to use. If present, \code{func} is
#' ignored.
#' @param WAIC Deprecated: If \code{TRUE}, use \code{func} to compute weights.
#' Otherwise tries to use DIC.
#' @param refresh Progress update refresh interval. 0 suppresses output.
#' @param replace Optional named list with replacement posterior samples. Used
#' for maginalizing over random effects, for example. See example in
#' \code{\link{link}}.
#' @param do_link If \code{TRUE}, compute and return \code{link} results
#' @param do_sim If \code{TRUE}, compute and return \code{sim} results
#' @author Richard McElreath
#' @seealso \code{\link{link}}, \code{\link{sim}}, \code{\link{compare}}
#' @export
ensemble <- function( ... , data , n=1e3 , func=WAIC , weights , refresh=0 , replace=list() , do_link=TRUE , do_sim=TRUE ) {
    # retrieve list of models
    L <- list(...)
    if ( is.list(L[[1]]) && length(L)==1 )
        L <- L[[1]]
    # retrieve model names from function call
    mnames <- match.call()
    mnames <- as.character(mnames)[2:(length(L)+1)]

    if ( missing(weights) ) {
        if ( length(L)>1 ) {
            use_func <- func
            ictab <- compare( ... , func=use_func , refresh=refresh , n=n , sort=FALSE )
            rownames(ictab@output) <- mnames
            weights <- ictab@output$weight
        } else {
            ictab <- NA
            weights <- 1
        }
    } else {
        # explicit custom weights
        ictab <- NA
        weights <- weights/sum(weights) # ensure sum to one
    }

    # compute number of predictions per model
    idx <- round( weights * n )

    # generate simulated predictions for each model
    # then compose ensemble using weights
    link.list <- list("empty")
    sim.list <- list("empty")
    if ( missing(data) ) {
        if ( do_link==TRUE )
            link.list <- lapply( L , function(m) link(m , n=n , refresh=refresh , replace=replace ) )
        if ( do_sim==TRUE )
            sim.list <- lapply( L , function(m) sim(m , n=n , refresh=refresh , replace=replace ) )
    } else {
        if ( do_link==TRUE )
            link.list <- lapply( L , function(m) link(m , data=data , n=n , refresh=refresh , replace=replace ) )
        if ( do_sim==TRUE )
            sim.list <- lapply( L , function(m) sim(m , data=data , n=n , refresh=refresh , replace=replace ) )
    }
    #print(str(sim.list))

    # compose weighted predictions
    # calculate indices corresponding to each model
    idx_start <- rep(1,length(idx))
    idx_end <- idx
    if ( length(L)>1 )
        for ( i in 2:length(idx) ) {
            idx_start[i] <- min( idx_start[i-1] + idx[i-1] , n )
            idx_end[i] <- min( idx_start[i] + idx[i] - 1 , n )
            if ( i==length(idx) ) idx_end[i] <- n
        }
    #print(idx)
    #print(idx_start)
    #print(idx_end)
    link_out <- link.list[[1]]
    sim_out <- sim.list[[1]]
    for ( i in 1:length(idx) ) {
        if ( idx[i]>0 ) {
            idxrange <- idx_start[i]:idx_end[i]
            if ( do_link==TRUE )
                link_out[idxrange,] <- link.list[[i]][idxrange,]
            if ( do_sim==TRUE )
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

