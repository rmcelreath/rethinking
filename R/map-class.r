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
    require(MASS)
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
    require(MASS)
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
