# DIC generic and methods

setGeneric("DIC",
function( object , ... ) {
    warning(paste("No specific DIC method for object of class ",class(object),". Returning AIC instead.\n",collapse="",sep=""))
    aic <- AIC(object)
    pD <- attr( logLik(object) , "df" )
    attr(aic,"pD") <- pD
    return(aic)
}
)

setMethod("DIC", "map",
function( object , n=1000 , ... ) {
    post <- extract.samples( object , n=n )
    dev <- sapply( 1:nrow(post) , 
        function(i) {
            p <- post[i,]
            names(p) <- names(post)
            2*object@fminuslogl( p ) 
        }
    )
    dev.hat <- deviance(object)
    val <- as.numeric( dev.hat + 2*( mean(dev) - dev.hat ) )
    attr(val,"pD") <- as.numeric( ( val - dev.hat )/2 )
    return(val)
}
)

setMethod("DIC", "map2stan",
function( object , ... ) {
    doDICmap2stan <- function(object) {
        fit <- object@stanfit
        # compute DIC
        dev.post <- extract(fit, "dev", permuted = TRUE, inc_warmup = FALSE)
        dbar <- mean( dev.post$dev )
        # to compute dhat, need to feed parameter averages back into compiled stan model
        post <- extract( fit )
        Epost <- list()
        for ( i in 1:length(post) ) {
            dims <- length( dim( post[[i]] ) )
            name <- names(post)[i]
            if ( name!="lp__" & name!="dev" ) {
                if ( dims==1 ) {
                    Epost[[ name ]] <- mean( post[[i]] )
                } else {
                    Epost[[ name ]] <- apply( post[[i]] , 2:dims , mean )
                }
            }
        }#i
        
        # push expected values back through model and fetch deviance
        #message("Taking one more sample now, at expected values of parameters, in order to compute DIC")
        fit2 <- stan( fit=fit , init=list(Epost) , data=object@data , pars="dev" , chains=1 , iter=1 , refresh=-1 )
        dhat <- as.numeric( extract(fit2,"dev") )
        pD <- dbar - dhat
        dic <- dbar + pD
        return( c( dic , pD ) )
    }
    
    # main
    if ( !is.null(attr(object,"DIC")) ) {
        val <- attr(object,"DIC")
        attr(val,"pD") <- attr(object,"pD")
    } else {
        # must compute
        v <- doDICmap2stan(object)
        val <- v[1]
        attr(val,"pD") <- v[2]
    }
    return(val)
}
)
