#' An S4 class to represent a MAP fit by quadratic approximation using \code{\link{map}()}
#'
#' @slot call The original call to \code{\link{map}()}
#' @slot coef The estimated posterior modes
#' @slot vcov Variance-covariance matrix
#' @slot optim List returned by \code{\link{optim}()}
#' @slot data The data
#' @slot formula Formula list
#' @slot formula_parsed Parsed version of \code{formula} as a list
#' of character.
#' @slot fminuslogl Negative log-likelihood function that accepts a vector of
#' parameter values
#' @slot links List of links
#' @slot start List of starting values for search algorithm
#'
#' @export
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

#' @export
setMethod("coef", "map", function(object) {
    object@coef
})

#' @export
setMethod("vcov", "map", function (object, ...) { object@vcov } )

#' @export
setMethod("nobs", "map", function (object, ...) { attr(object,"nobs") } )

#' @export
setGeneric("stancode",
function( object , ... ) {
    print(class(object))
}
)
#' @export
setMethod("stancode", "map",
function(object) {
    # compile Stan code through map2stan
    m <- map2stan( object@formula , data=object@data , start=object@start , sample=FALSE )
    # display
    cat( m$model )
    return( invisible( m$model ) )
}
)

#' @export
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

#' @export
setMethod("deviance", "map",
function (object, ...)
{
  -2*logLik(object)
})

#' @export
setMethod("AIC", "map",
function (object, ...)
{
  -2*logLik(object) + 2*length(coef(object))
})

#' @export
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

#' @export
setMethod("summary", "map", function(object){

    precis(object)

})

#' Collect posterior samples from a fit model
#'
#' Extracts or draws posterior samples for fit models.
#'
#' For \code{map2stan} models, this function returns cleaned samples already
#' contained in the object. These samples are cleaned of dimension attributes
#' and the \code{lp__} and \code{dev} traces that are used internally.
#'
#' For \code{map} and other types, it uses the variance-covariance matrix and
#' coefficients to define a multivariate Gaussian posterior to draw \code{n}
#' samples from.
#'
#' @param object Fit model to extract samples from
#' @param n Number of samples to simulate
#' @param clean.names When \code{TRUE}, removes parentheses from parameters
#' named \code{(Intercept)}.
#' @param ... Other parameters to pass to descendent functions (when defined)
#' @return A named \code{list} (for \code{map2stan}) or \code{data.frame}
#' containing samples for each parameter in the posterior distribution.
#' @author Richard McElreath
#' @seealso \code{\link{mvrnorm}}

#' @export

setGeneric("extract.samples",
function( object , n=10000 , clean.names=TRUE , ... ) {
    #require(MASS)
    mu <- 0
    if ( class(object)[1] %in% c("mer","bmer","glmerMod","lmerMod") ) {
        mu <- fixef(object)
    } else {
        mu <- xcoef(object)
    }
    result <- as.data.frame( MASS::mvrnorm( n=n , mu=mu , Sigma=vcov(object) ) )
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

#' @export
setMethod("extract.samples", "map",
function(object,n=1e4,...){
    # require(MASS) # import now, so no need to require?
    mu <- object@coef
    result <- as.data.frame( MASS::mvrnorm( n=n , mu=mu , Sigma=vcov(object) ) )
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


#' @export
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

