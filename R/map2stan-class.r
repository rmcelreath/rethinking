setClass("map2stan", representation( call = "language",
                                stanfit = "stanfit",
                                coef = "numeric",
                                vcov = "matrix",
                                data = "list",
                                start = "list",
                                formula = "list" ,
                                formula_parsed = "list" ))

setMethod("coef", "map2stan", function(object) {
    object@coef
})

extract.samples <- function(object) {
    p <- extract(object@stanfit)
    # get rid of dev and lp__
    p[['dev']] <- NULL
    p[['lp__']] <- NULL
    # get rid of those ugly dimnames
    for ( i in 1:length(p) ) {
        attr(p[[i]],"dimnames") <- NULL
    }
    return(p)
}

plotchains <- function(object,...) {
    if ( class(object)=="map2stan" )
        rstan::traceplot( object@stanfit , ask=TRUE , pars=names(object@start) , ... )
}

plotpost <- function(object,n=1000,col=col.alpha("slateblue",0.3),cex=0.8,pch=16,...) {
    o <- as.data.frame(object)
    pairs(o[1:n,],col=col,cex=cex,pch=pch,...)
}

stancode <- function(object) {
    cat( object@stanfit@stanmodel@model_code )
    return( invisible( object@stanfit@stanmodel@model_code ) )
}

setMethod("vcov", "map2stan", function (object, ...) { object@vcov } )

setMethod("nobs", "map2stan", function (object, ...) { attr(object,"nobs") } )

setMethod("logLik", "map2stan",
function (object, ...)
{
    if(length(list(...)))
        warning("extra arguments discarded")
    val <- (-1)*attr(object,"deviance")/2
    attr(val, "df") <- length(object@coef)
    attr(val, "nobs") <- attr(object,"nobs")
    class(val) <- "logLik"
    val
  })
  
setMethod("deviance", "map2stan",
function (object, ...)
{
  attr(object,"deviance")
})

DIC <- function( object ) {
    val <- attr(object,"DIC")
    attr(val,"pD") <- attr(object,"pD")
    val
}

setMethod("show", "map2stan", function(object){

    cat("\nmap2stan model fit\n")
    
    cat("\nFormula:\n")
    for ( i in 1:length(object@formula) ) {
        print( object@formula[[i]] )
    }
    
    #cat("\nExpected values of fixed effects:\n")
    #print(coef(object))
    
    cat("\nLog-likelihood at expected values: ")
    cat(round(as.numeric(logLik(object)),2),"\n")
    
    cat("Deviance: ")
    cat(round(as.numeric(deviance(object)),2),"\n")
    
    cat("DIC: ")
    cat(round(as.numeric(DIC(object)),2),"\n")
    
    cat("Effective number of parameters (pD): ")
    cat(round(as.numeric(attr(object,"pD")),2),"\n")
    
  })

setMethod("summary", "map2stan", function(object){
    
    show(object@stanfit)
    
})