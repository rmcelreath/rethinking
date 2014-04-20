setClass("map", representation( call = "language",
                                coef = "numeric",
                                vcov = "matrix",
                                optim = "list",
                                data = "list",
                                formula = "list",
                                formula_parsed = "list",
                                fminuslogl = "function",
                                links = "list"))

setMethod("coef", "map", function(object) {
    object@coef
})

setMethod("vcov", "map", function (object, ...) { object@vcov } )

setMethod("nobs", "map", function (object, ...) { attr(object,"nobs") } )

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