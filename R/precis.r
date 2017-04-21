# my model summary function, précis

precis.whitelist <- data.frame(stringsAsFactors = FALSE,
    class=c("map","map2stan","lm","glm","mle2","mer","bmer","polr","data.frame","clmm","clmm2","list","stanfit","lmerMod","glmerMod") ,
    coef.method=c("coef","coef","coef","coef","coef","fixef.plus","fixef.plus","polr","chain","coef","coef","mcarray","stanfit","fixef.plus","fixef.plus") ,
    vcov.method=c("vcov","vcov","vcov","vcov","vcov","vcov.VarCorr","vcov.VarCorr","vcov","chain","vcov","vcov","mcarray","stanfit","vcov.VarCorr","vcov.VarCorr") ,
    nobs.method=c("nobs","nobs","nobs","nobs","mle2","mer","mer","nobs","chain","nobs","nobs","chain","stanfit","mer","mer")
)

# precis class definition and show method
#' @export
setClass( "precis" , representation( output="data.frame" , digits="numeric" ) )

precis_show <- function( object ) {
    #print( round( object@output , object@digits ) )
    r <- format_show( object@output , digits=c('default__'=object@digits,'n_eff'=0) )
    print(r)
}

#' @export
setMethod( "show" , "precis" , function(object) precis_show(object) )

precis_plot <- function( x , y , pars , col.ci="black" , xlab="Value" , ... ) {
    x <- x@output
    if ( !missing(pars) ) {
        x <- x[pars,]
    }
    n <- nrow(x)
    mu <- x[n:1,1]
    left <- x[[3]][n:1]
    right <- x[[4]][n:1]
    set_nice_margins()
    dotchart( mu , labels=rownames(x)[n:1] , xlab=xlab , xlim=c(min(left),max(right)) , ... )
    for ( i in 1:length(mu) ) lines( c(left[i],right[i]) , c(i,i) , lwd=2 , col=col.ci )
    abline( v=0 , lty=1 , col=col.alpha("black",0.15) )
}

#' ggplot2 plot of precis of models fit with rethinking package
#'
#' @param x a precis object
#' @param params character vectro of parameters to include
#' @param col_ci color of interval
#' @param dot_size size of dots
#' @param line_size size of lines
#' @param flip logical indicating whther x- and y- axes should be swapped.
#'   If TRUE (the default), then parameters will be listed along the y-axis.
#' @param ... additional arguments (currently ignored)
#' @export
#'
ggplot.precis <- function( x, params, col_ci = "black", lab = "value", flip = TRUE,
                           dot_size = 1, line_size = 0.5, ...) {
  x <- x@output
  if ( !missing(params) ) {
    x <- x[params,]
  }
  n <- nrow(x)
  PlotData <- data_frame(
    mu = x[n:1,1],
    left = x[[3]][n:1],
    right = x[[4]][n:1],
    param_name = factor(rownames(x)[n:1], labels = rownames(x)[n:1])
  )
  p <-
    ggplot(PlotData) +
    geom_hline(yintercept = 0, color = "red", alpha = 0.2) +
    geom_linerange(
      aes(y = mu, x = param_name, ymin = left, ymax = right),
      size = line_size) +
    geom_point(aes(y = mu, x = param_name), size = dot_size) +
    labs(y = lab, x = "parameter")

  if (flip) p <- p + coord_flip()
  p
}
#' @export
setMethod( "plot" , "precis" , function(x,y,...) precis_plot(x,y,...) )

# function to process a list of posterior samples from extract.samples into a summary table
# needed because as.data.frame borks the ordering of matrix parameters like varying effects
postlistprecis <- function( post , prob=0.95 ) {
    n_pars <- length(post)
    result <- data.frame( Mean=0 , StdDev=0 , lower=0 , upper=0 )
    r <- 1
    for ( k in 1:n_pars ) {
        dims <- dim( post[[k]] )
        if ( length(dims)==1 ) {
            # single parameter
            hpd <- as.numeric( HPDI( post[[k]] , prob=prob ) )
            result[r,] <- c( mean(post[[k]]) , sd(post[[k]]) , hpd[1] , hpd[2] )
            rownames(result)[r] <- names(post)[k]
            r <- r + 1
        }
        if ( length(dims)==2 ) {
            # vector of parameters
            # loop over
            for ( i in 1:dims[2] ) {
                hpd <- as.numeric( HPDI( post[[k]][,i] , prob=prob ) )
                result[r,] <- c( mean(post[[k]][,i]) , sd(post[[k]][,i]) , hpd[1] , hpd[2] )
                rownames(result)[r] <- concat( names(post)[k] , "[" , i , "]" )
                r <- r + 1
            }
        }
        if ( length(dims)==3 ) {
            # matrix of parameters
            for ( i in 1:dims[2] ) {
                for ( j in 1:dims[3] ) {
                    hpd <- as.numeric( HPDI( post[[k]][,i,j] , prob=prob ) )
                    result[r,] <- c( mean(post[[k]][,i,j]) , sd(post[[k]][,i,j]) , hpd[1] , hpd[2] )
                    rownames(result)[r] <- concat( names(post)[k] , "[" , i , "," , j , "]" )
                    r <- r + 1
                }
            }
        }
    }
    colnames(result)[3:4] <- c(paste("lower", prob), paste("upper", prob))
    result
}



#' Precis of model fit
#'
#' Displays concise parameter estimate information for an existing model fit.
#'
#' Creates a table of estimates and standard errors, with optional confidence
#' intervals and parameter correlations. Posterior intervals are quadratic
#' estimates, derived from standard deviations, unless the model uses samples
#' from the posterior distribution, in which case \code{\link{HPDI}} is used
#' instead.
#'
#' Can also provide expected value, standard deviation, and HPDI columns for a
#' data frame.
#'
#' @param model Fit model object
#' @param depth If \code{1}, suppresses vectors and matrices of parameters. If
#' \code{2}, displays all parameters
#' @param pars Optional character vector of parameter names to display
#' @param ci Show quadratic estimate confidence intervals
#' @param prob Width of posterior intervals
#' @param corr If \code{TRUE}, show correlations among parameters in output
#' @param digits Number of decimal places to display in output
#' @param warn If \code{TRUE}, warns about various things
#' @return A data frame with a row for each parameter.
#' @author Richard McElreath
#' @export
precis <- function( model , depth=1 , pars , ci=TRUE , prob=0.89 , corr=FALSE , digits=2 , warn=TRUE ) {
    the.class <- class(model)[1]
    found.class <- FALSE
    if ( inherits(model, "numeric") ) {
        # single vector of values
        # coerce to data frame
        model <- as.data.frame(model)
        the.class <- class(model)[1]
    }
    if ( inherits(model, precis.whitelist$class) ) found.class <- TRUE
    if ( inherits(model, "list") )
        if ( ! inherits(model[[1]], "mcarray") ) found.class <- FALSE
    if ( found.class ) {
        est <- xcoef( model )
        se <- xse( model )
        if ( corr==TRUE ) Rho <- xrho( model )
    }
    if ( ! found.class ) {
        message("No handler found for model (bad class)")
        return(invisible())
    }
    # format
    fname <- "Mean"
    # capitalize first letter
    fname <- concat( toupper(substring(fname,1,1)) , substring(fname,2) )
    result <- data.frame( est=est , se=se )
    colnames(result) <- c( fname ,"StdDev")
    if ( ci ) {
        ci <- confint_quad( est=est , se=se , prob=prob )
        if ( inherits(model, "data.frame") ) {
            # HPDI from samples
            ci <- t( apply( model , 2 , HPDI , prob=prob ) )
        }
        result <- cbind( result , ci )
        if ( inherits(model, "map2stan") ) {
            # HPDI from samples
            post <- extract.samples(model)
            result <- postlistprecis( post , prob=prob )
        }
        if ( inherits(model, "stanfit") ) {
            # HPDI from samples
            post <- extract.samples(model)
            post[['lp__']] <- NULL
            result <- postlistprecis( post , prob=prob )
        }
    }
    if ( inherits(model, c("map2stan", "stanfit")) ) {
        # add n_eff to result
        #require(rstan)
        if ( inherits(model, "map2stan") )
            the_summary <- summary( model@stanfit )$summary
        else
            the_summary <- summary( model )$summary
        n_eff <- the_summary[,'n_eff']
        n_eff <- n_eff[ -which(names(n_eff)=="lp__") ]
        Rhat <- the_summary[,'Rhat']
        Rhat <- Rhat[ -which(names(Rhat)=="lp__") ]
        if ( inherits(model, "map2stan") ) {
            n_eff <- n_eff[ -which(names(n_eff)=="dev") ]
            Rhat <- Rhat[ -which(names(Rhat)=="dev") ]
        }
        result <- cbind( result , n_eff , Rhat )

        # check divergent iterations
        nd <- divergent(model)
        if ( nd > 0 & warn==TRUE ) {
            warning( concat("There were ",nd," divergent iterations during sampling.\nCheck the chains (trace plots, n_eff, Rhat) carefully to ensure they are valid.\n") )
        }
    }
    if ( corr ) {
        result <- cbind( result , Rho )
    }

    #if ( type.s )
    #    result[,"Pr(S)"] <- format.pval( type.s( est , se ) )

    if ( precis.whitelist$vcov.method[ precis.whitelist$class==the.class ]=="vcov.VarCorr" ) {
        message( "Quadratic approximation (standard errors) unreliable for variance components. Use MCMC to estimate precision of variance components." )
    }

    # deal with depth
    if ( depth == 1 ) {
        hits <- regexpr("]",rownames(result),fixed=TRUE)
        hits_idx <- which( hits > -1 )
        if ( length(hits_idx)>0 ) {
            result <- result[-hits_idx,]
            message( paste( length(hits_idx) , "vector or matrix parameters omitted in display. Use depth=2 to show them." ) )
        }
    }

    # deal with pars list
    if ( !missing(pars) ) {
        # have to handle vector/matrix parameters,
        # so split off any [.] and then prune with names only
        clean_names <- as.character( sapply( rownames(result) , function(s) strsplit( s , "[" , fixed=TRUE )[[1]][1] ) )
        result <- result[ clean_names %in% pars , ]
    }

    # result
    new( "precis" , output=result , digits=digits )
}

####
xcoef <- function( model , func=mean ) {
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
        result <- apply( model , 2 , func )
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
        if ( !is.na(sigma.resid) ) result['(residual)'] <- sigma.resid
    }
    xcheckconvergence( model )
    result
}

xcheckconvergence <- function( model ) {
    the.class <- class(model)[1]
    k <- 0
    if ( the.class=="mle2" ) {
        if ( model@details$convergence != 0 ) {
            k <- model@details$convergence
        }
    }
    if ( the.class=="map" ) {
        if ( model@optim$convergence != 0 ) {
            k <- model@optim$convergence
        }
    }

    if ( k > 0 ) {
        message( paste("Caution, model may not have converged.") )
        if ( k==1 ) {
            message( "Code 1: Maximum iterations reached." )
        }
        if ( k==10 ) {
            message( "Code 10: Degenerate Nelder-Mead simplex." )
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
#' @export
rrformat <- function( matrix , digits=2 , width=7 ) {
    if ( length(digits)==1 ) digits <- rep(digits,nrow(matrix))
    result <- matrix
    for ( i in 1:nrow(matrix) ) {
        result[i,] <- format( round(matrix[i,],digits[i]) , width=width )
    }
    result
}

