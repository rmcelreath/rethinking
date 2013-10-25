# my model summary function, précis

precis.whitelist <- data.frame( 
    class=c("map","lm","glm","mle2","mer","bmer","polr","data.frame","clmm","clmm2","list","stanfit","lmerMod","glmerMod") , 
    coef.method=c("coef","coef","coef","coef","fixef.plus","fixef.plus","polr","chain","coef","coef","mcarray","stanfit","fixef.plus","fixef.plus") , 
    vcov.method=c("vcov","vcov","vcov","vcov","vcov.VarCorr","vcov.VarCorr","vcov","chain","vcov","vcov","mcarray","stanfit","vcov.VarCorr","vcov.VarCorr") ,
    nobs.method=c("nobs","nobs","nobs","mle2","mer","mer","nobs","chain","nobs","nobs","chain","stanfit","mer","mer")
)

# precis class definition and show method
setClass( "precis" , representation( output="data.frame" , digits="numeric" ) )
precis.show <- function( object ) {
    print( round( object@output , object@digits ) )
}
setMethod( "show" , "precis" , function(object) precis.show(object) )

precis.plot <- function( x , y , pars , col.ci="black" , ... ) {
    x <- x@output
    if ( !missing(pars) ) {
        x <- x[pars,]
    }
    n <- nrow(x)
    mu <- x$Estimate[n:1]
    left <- x[[3]][n:1]
    right <- x[[4]][n:1]
    dotchart( mu , labels=rownames(x)[n:1] , xlab="Estimate" , xlim=c(min(left),max(right)) , ... )
    for ( i in 1:length(mu) ) lines( c(left[i],right[i]) , c(i,i) , lwd=2 , col=col.ci )
    abline( v=0 , lty=1 , col=col.alpha("black",0.15) )
}
setMethod( "plot" , "precis" , function(x,y,...) precis.plot(x,y,...) )

precis <- function( model , type.s=FALSE , ci=TRUE , level=0.95 , corr=FALSE , digits=2 , warn=TRUE ) {
    the.class <- class(model)[1]
    found.class <- FALSE
    if ( the.class=="numeric" ) {
        # single vector of values
        # coerce to data frame
        model <- as.data.frame(model)
        the.class <- class(model)[1]
    }
    if ( any( precis.whitelist$class==the.class ) ) found.class <- TRUE
    if ( the.class=="list" )
        if ( class( model[[1]] ) != "mcarray" ) found.class <- FALSE
    if ( found.class==TRUE ) {
        est <- xcoef( model )
        se <- xse( model )
        if ( corr==TRUE ) Rho <- xrho( model )
    }
    if ( found.class==FALSE ) {
        return( paste("No handler found for model of class",the.class) )
    }
    # format
    result <- data.frame( est=est , se=se )
    colnames(result) <- c("Estimate","S.E.")
    if ( ci==TRUE ) {
        ci <- confint.quad( est=est , se=se , level=level )
        if ( the.class=="data.frame" ) {
            # PCI from chain
            ci <- t( apply( model , 2 , PCI , prob=level ) )
            colnames(result)[1] <- "Expectation"
            colnames(result)[2] <- "StdDev"
        }
        result <- cbind( result , ci )
    }
    if ( corr==TRUE ) {
        result <- cbind( result , Rho )
    }
    if ( type.s==TRUE )
        result[,"Pr(S)"] <- format.pval( type.s( est , se ) )
    if ( precis.whitelist$vcov.method[ precis.whitelist$class==the.class ]=="vcov.VarCorr" ) {
        message( "Quadratic approximation (standard errors) unreliable for variance components. Use MCMC to estimate precision of variance components." )
    }
    new( "precis" , output=result , digits=digits )
}

####
xcoef <- function( model ) {
    the.class <- class(model)[1]
    the.method <- precis.whitelist$coef.method[ precis.whitelist$class==the.class ]
    if ( the.method=="coef" ) {
        result <- coef(model)
    }
    if ( the.method=="fixef" ) {
        result <- fixef(model)
    }
    if ( the.method=="polr" ) {
        result <- summary(model)$coefficients[,1]
    }
    if ( the.method=="chain" ) {
        # average of chains
        result <- apply( model , 2 , mean )
    }
    if ( the.method=="stanfit" ) {
        result <- summary( model )$summary[,1]
    }
    if ( the.method=="mcarray" ) {
        # jags.samples result, hopefully
        result <- NULL
        result.names <- NULL
        for ( j in 1:length(model) ) {
            # explode compact arrays of coefficients
            dims <- dim( model[[j]] )
            if ( length(dims)==3 ) {
                est <- rep(0,dims[1])
                for ( k in 1:dims[1] ) {
                    est[k] <- mean( as.vector(model[[j]][ k , , ]) ) # marginalize over iteration and chain
                }
                result <- c( result , est )
                if ( dims[1] > 1 ) {
                    newnames <- paste( names(model)[j] , "[" , 1:dims[1] , "]" , sep="" )
                } else {
                    newnames <- names(model)[j]
                }
                result.names <- c( result.names , newnames )
            } # if dim 3
        } # for each array
        names(result) <- result.names
    }
    if ( the.method=="fixef.plus" ) {
        # fixef from lmer plus variance components
        result <- fixef(model)
        vc <- VarCorr(model)
        clusters <- names(vc)
        for( i in 1:length(clusters) ) {
            sigma <- sqrt(diag(vc[[i]]))
            names(sigma) <- paste( rownames(vc[[i]]) , clusters[i] , sep="|" )
            names(sigma) <- paste( "(" , names(sigma) , ")" , sep="" )
            result <- c( result , sigma )
        }
        sigma.resid <- attr( vc , "sc" )
        # new lme4 uses "useSc" attribute to flag for residual variance
        # so "sc" is always present, but "useSc" may be TRUE or FALSE
        # should worry about deprecated lme4 model fits?
        if ( !is.na(sigma.resid) & attr( vc , "useSc" ) ) {
            names(sigma.resid) <- "(Residual)"
            result <- c( result , sigma.resid )
        }
    }
    xcheckconvergence( model )
    result
}

xcheckconvergence <- function( model ) {
    the.class <- class(model)[1]
    if ( the.class=="mle2" ) {
        if ( model@details$convergence != 0 ) {
            k <- model@details$convergence
            message( paste("Caution, model may not have converged.") )
            if ( k==1 ) {
                message( "Code 1: Maximum iterations reached." )
            }
            if ( k==10 ) {
                message( "Code 10: Degenerate Nelder-Mead simplex." )
            }
        }
    }
}

xse <- function( model ) {
    the.class <- class(model)[1]
    the.method <- precis.whitelist$vcov.method[ precis.whitelist$class==the.class ]
    if ( the.method=="vcov" ) {
        result <- sqrt(diag(vcov(model)))
    }
    if ( the.method=="vcov.VarCorr" ) {
        result <- sqrt(diag( as.matrix(vcov(model)) ))
        num.sigma <- length( xcoef(model) ) - length( fixef(model) )
        result <- c( result , rep(NA,num.sigma) )
    }
    if ( the.method=="chain" ) {
        # sd of chains
        result <- apply( model , 2 , sd )
    }
    if ( the.method=="stanfit" ) {
        result <- summary( model )$summary[,3]
    }
    if ( the.method=="mcarray" ) {
        # jags.samples result, hopefully
        result <- NULL
        result.names <- NULL
        for ( j in 1:length(model) ) {
            # explode compact arrays of coefficients
            dims <- dim( model[[j]] )
            if ( length(dims)==3 ) {
                est <- rep(0,dims[1])
                for ( k in 1:dims[1] ) {
                    est[k] <- sd( as.vector(model[[j]][ k , , ]) ) # marginalize over iteration and chain
                }
                result <- c( result , est )
                if ( dims[1] > 1 ) {
                    newnames <- paste( names(model)[j] , "[" , 1:dims[1] , "]" , sep="" )
                } else {
                    newnames <- names(model)[j]
                }
                result.names <- c( result.names , newnames )
            } # if dim 3
        } # for each array
        names(result) <- result.names
    }
    result
}

xrho <- function( model ) {
    the.class <- class(model)[1]
    the.method <- precis.whitelist$vcov.method[ precis.whitelist$class==the.class ]
    if ( the.method=="vcov" ) {
        result <- cov2cor( vcov(model) ) 
    }
    if ( the.method=="vcov.VarCorr" ) {
        result <- sqrt(diag( as.matrix(vcov(model)) ))
        num.sigma <- length( xcoef(model) ) - length( fixef(model) )
        result <- c( result , rep(NA,num.sigma) )
    }
    if ( the.method=="chain" ) {
        # sd of chains
        result <- apply( model , 2 , sd )
    }
    if ( the.method=="stanfit" ) {
        result <- summary( model )$summary[,3]
    }
    if ( the.method=="mcarray" ) {
        # jags.samples result, hopefully
        result <- NULL
        result.names <- NULL
        for ( j in 1:length(model) ) {
            # explode compact arrays of coefficients
            dims <- dim( model[[j]] )
            if ( length(dims)==3 ) {
                est <- rep(0,dims[1])
                for ( k in 1:dims[1] ) {
                    est[k] <- sd( as.vector(model[[j]][ k , , ]) ) # marginalize over iteration and chain
                }
                result <- c( result , est )
                if ( dims[1] > 1 ) {
                    newnames <- paste( names(model)[j] , "[" , 1:dims[1] , "]" , sep="" )
                } else {
                    newnames <- names(model)[j]
                }
                result.names <- c( result.names , newnames )
            } # if dim 3
        } # for each array
        names(result) <- result.names
    }
    result
}

xnobs <- function( model ) {
    the.class <- class(model)[1]
    the.method <- precis.whitelist$nobs.method[ precis.whitelist$class==the.class ]
    if ( the.method=="nobs" ) result <- nobs(model)
    if ( the.method=="mle2" ) result <- length(model@data[[1]])
    if ( the.method=="mer" ) result <- nrow(model@frame)
    if ( the.method=="chain" ) result <- nrow(model)
    if ( the.method=="stanfit" ) result <- 0
    result
}

# row-by-row matrix formatting function
rrformat <- function( matrix , digits=2 , width=7 ) {
    if ( length(digits)==1 ) digits <- rep(digits,nrow(matrix))
    result <- matrix
    for ( i in 1:nrow(matrix) ) {
        result[i,] <- format( round(matrix[i,],digits[i]) , width=width )
    }
    result
}

