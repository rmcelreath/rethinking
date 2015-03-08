setClass("map", representation( call = "language",
                                coef = "numeric",
                                vcov = "matrix",
                                optim = "list",
                                data = "list",
                                start = "list",
                                formula = "list",
                                formula_parsed = "list",
                                fminuslogl = "function",
                                links = "list"))

setMethod("coef", "map", function(object) {
    object@coef
})

setMethod("vcov", "map", function (object, ...) { object@vcov } )

setMethod("nobs", "map", function (object, ...) { attr(object,"nobs") } )

setGeneric("stancode",
function( object , ... ) {
    print(class(object))
}
)
setMethod("stancode", "map",
function(object) {
    # compile Stan code through map2stan
    m <- map2stan( object@formula , data=object@data , start=object@start , sample=FALSE )
    # display
    cat( m$model )
    return( invisible( m$model ) )
}
)

setMethod("logLik", "map",
function (object, ...)
{
    if(length(list(...)))
        warning("extra arguments discarded")
    val <- -object@optim[['minuslogl']]
    attr(val, "df") <- length(object@coef)
    attr(val, "nobs") <- attr(object,"nobs")
    class(val) <- "logLik"
    val
  })
  
setMethod("deviance", "map",
function (object, ...)
{
  -2*logLik(object)
})

setMethod("AIC", "map",
function (object, ...)
{
  -2*logLik(object) + 2*length(coef(object))
})

setMethod("show", "map", function(object){

    cat("\nMaximum a posteriori (MAP) model fit\n")
    
    cat("\nFormula:\n")
    for ( i in 1:length(object@formula) ) {
        print( object@formula[[i]] )
    }
    
    cat("\nMAP values:\n")
    print(coef(object))
    
    cat("\nLog-likelihood: ")
    cat(round(as.numeric(logLik(object)),2),"\n")
    
    if ( object@optim$convergence > 0 )
      cat("\nWarning: optimization did not converge (code ",
          object@optim$convergence, ": " , object@optim$message,")\n",sep="")
  })
  
setMethod("summary", "map", function(object){
    
    precis(object)
    
})

setGeneric("extract.samples",
function( object , n=10000 , clean.names=TRUE , ... ) {
    #require(MASS)
    mu <- 0
    if ( class(object)[1] %in% c("mer","bmer","glmerMod","lmerMod") ) {
        mu <- fixef(object)
    } else {
        mu <- xcoef(object)
    }
    result <- as.data.frame( mvrnorm( n=n , mu=mu , Sigma=vcov(object) ) )
    if ( clean.names==TRUE ) {
        # convert (Intercept) to Intercept
        for ( i in 1:ncol(result) ) {
            if ( colnames(result)[i] == "(Intercept)" ) {
                colnames(result)[i] <- "Intercept"
            }
        }
    }
    result
}
)

setMethod("extract.samples", "map",
function(object,n=1e4,...){
    # require(MASS) # import now, so no need to require?
    mu <- object@coef
    result <- as.data.frame( mvrnorm( n=n , mu=mu , Sigma=vcov(object) ) )
    # convert vector parameters to vectors in list
    veclist <- attr(object,"veclist")
    name_head <- function(aname) strsplit( aname , "[" , fixed=TRUE )[[1]][1]
    name_index <- function(aname) as.numeric(regmatches( aname , regexec( "\\[(.+)\\]" , aname ) )[[1]][2])
    if ( length(veclist) > 0 ) {
        new_result <- list()
        # copy non-vector samples into new list
        for ( i in 1:length(result) ) {
            if ( !( name_head(names(result)[i]) %in% names(veclist) ) ) {
                new_result[[ names(result)[i] ]] <- result[[i]]
            }
        }#i
        # now build vectors out of paramters with [n] in name
        for ( i in 1:length(veclist) ) {
            # empty matrix with parameters on cols and samples on rows
            # so n-by-m, where m in number of pars in vector and n is number of samples
            new_matrix <- matrix( 0 , ncol=veclist[[i]]$n , nrow=n )
            for ( j in 1:length(result) ) {
                if ( name_head(names(result)[j]) == names(veclist)[i] ) {
                    the_index <- name_index( names(result)[j] )
                    new_matrix[,the_index] <- result[[j]]
                }
            }
            new_result[[ names(veclist)[i] ]] <- new_matrix
        }#i
        result <- new_result
    }
    # return result
    result
})


setMethod("pairs" , "map" , function(x, n=500 , alpha=0.7 , cex=0.7 , pch=16 , adj=1 , pars , ...) {
    posterior <- extract.samples(x)
    if ( !missing(pars) ) {
        # select out named parameters
        p <- list()
        for ( k in pars ) p[[k]] <- posterior[[k]]
        posterior <- p
    }
    panel.dens <- function(x, ...) {
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(usr[1:2], 0, 1.5) )
        h <- density(x,adj=adj)
        y <- h$y
        y <- y/max(y)
        abline( v=0 , col="gray" , lwd=0.5 )
        lines( h$x , y )
    }
    panel.2d <- function( x , y , ... ) {
        i <- sample( 1:length(x) , size=n )
        abline( v=0 , col="gray" , lwd=0.5 )
        abline( h=0 , col="gray" , lwd=0.5 )
        dcols <- densCols( x[i] , y[i] )
        dcols <- sapply( dcols , function(k) col.alpha(k,alpha) )
        points( x[i] , y[i] , col=dcols , ... )
    }
    panel.cor <- function( x , y , ... ) {
        k <- cor( x , y )
        cx <- sum(range(x))/2
        cy <- sum(range(y))/2
        text( cx , cy , round(k,2) , cex=2*exp(abs(k))/exp(1) )
    }
    pairs( posterior , cex=cex , pch=pch , upper.panel=panel.2d , lower.panel=panel.cor , diag.panel=panel.dens , ... )
})

