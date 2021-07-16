# my model summary function, précis

precis.whitelist <- data.frame( 
    class=c("map","map2stan","lm","glm","mle2","mer","bmer","polr","data.frame","clmm","clmm2","list","stanfit","lmerMod","glmerMod","ulam") , 
    coef.method=c("coef","coef","coef","coef","coef","fixef.plus","fixef.plus","polr","chain","coef","coef","mcarray","stanfit","fixef.plus","fixef.plus","coef") , 
    vcov.method=c("vcov","vcov","vcov","vcov","vcov","vcov.VarCorr","vcov.VarCorr","vcov","chain","vcov","vcov","mcarray","stanfit","vcov.VarCorr","vcov.VarCorr","vcov") ,
    nobs.method=c("nobs","nobs","nobs","nobs","mle2","mer","mer","nobs","chain","nobs","nobs","chain","stanfit","mer","mer","nobs")
)

# precis class definition and show method
# setClass( "precis" , representation( output="data.frame" , digits="numeric" ) )
setClass( "precis" , slots=c( digits="numeric" ) , contains="data.frame" )

precis_show_old <- function( object ) {
    #print( round( object@output , object@digits ) )
    r <- format_show( object@output , digits=c('default__'=object@digits,'n_eff'=0) )
    print(r)
}
precis_show <- function( object ) {
    #print( round( object@output , object@digits ) )
    r <- format_show( object , digits=c('default__'=object@digits,'n_eff'=0) )
    has_header <- !is.null( attr(object,"header") )
    if ( has_header ) {
        # show header
        cat( attr(object,"header") )
        cat( "\n" )
    }
    print(r)
}

setMethod( "show" , "precis" , function(object) precis_show(object) )

precis_plot <- function( x , y , pars , col.ci="black" , xlab="Value" , add=FALSE , xlim=NULL , labels=rownames(x)[1:n] , ... ) {
    if ( !missing(pars) ) {
        x <- x[pars,]
    }
    n <- nrow(x)
    mu <- x[n:1,1]
    left <- x[[3]][n:1]
    right <- x[[4]][n:1]
    set_nice_margins()
    labels <- labels[n:1]
    if ( is.null(xlim) ) xlim <- c(min(left),max(right))
    if ( add==FALSE )
        dotchart( mu , labels=labels , xlab=xlab , xlim=xlim , ... )
    else
        points( mu[n:1] , n:1 , ... )
    for ( i in 1:length(mu) ) lines( c(left[i],right[i]) , c(i,i) , lwd=2 , col=col.ci )
    if ( add==FALSE ) abline( v=0 , lty=1 , col=col.alpha("black",0.15) )
}
setMethod( "plot" , "precis" , function(x,y,...) precis_plot(x,y,...) )

# function to process a list of posterior samples from extract.samples into a summary table
# needed because as.data.frame borks the ordering of matrix parameters like varying effects
postlistprecis <- function( post , prob=0.95 , spark=FALSE ) {
    n_pars <- length(post)
    result <- data.frame( Mean=0 , StdDev=0 , lower=0 , upper=0 )
    if ( spark!=FALSE ) {
        result <- data.frame( Mean=0 , StdDev=0 , Min=0 , Distribution="" , Max=0 )
        result$Distribution <- as.character(result$Distribution)
    }
    r <- 1
    for ( k in 1:n_pars ) {
        dims <- dim( post[[k]] )
        if ( length(dims)==1 ) {
            # single parameter
            if ( spark==FALSE ) {
                hpd <- as.numeric( HPDI( post[[k]] , prob=prob ) )
                result[r,] <- c( mean(post[[k]]) , sd(post[[k]]) , hpd[1] , hpd[2] )
            } else {
                # histosparks in place of HPDI
                the_spark <- histospark( post[[k]] , width=spark )
                result[r,1:3] <- c( mean(post[[k]]) , sd(post[[k]]) , min(post[[k]]) )
                result[r,4] <- the_spark
                result[r,5] <- max( post[[k]] )
            }
            rownames(result)[r] <- names(post)[k]
            r <- r + 1
        }
        if ( length(dims)==2 ) {
            # vector of parameters
            # loop over
            for ( i in 1:dims[2] ) {
                if ( spark==FALSE ) {
                    hpd <- as.numeric( HPDI( post[[k]][,i] , prob=prob ) )
                    result[r,] <- c( mean(post[[k]][,i]) , sd(post[[k]][,i]) , hpd[1] , hpd[2] )
                } else {
                    the_spark <- histospark( post[[k]][,i] , width=spark )
                    result[r,1:3] <- c( mean(post[[k]][,i]) , sd(post[[k]][,i]) , min(post[[k]][,i]) )
                    result[r,4] <- the_spark
                    result[r,5] <- max( post[[k]][,i] )
                }
                rownames(result)[r] <- concat( names(post)[k] , "[" , i , "]" )
                r <- r + 1
            }
        }
        if ( length(dims)==3 ) {
            # matrix of parameters
            for ( i in 1:dims[2] ) {
                for ( j in 1:dims[3] ) {
                    if ( spark==FALSE ) {
                        hpd <- as.numeric( HPDI( post[[k]][,i,j] , prob=prob ) )
                        result[r,] <- c( mean(post[[k]][,i,j]) , sd(post[[k]][,i,j]) , hpd[1] , hpd[2] )
                    } else {
                        the_spark <- histospark( post[[k]][,i,j] , width=spark )
                        result[r,1:3] <- c( mean(post[[k]][,i,j]) , sd(post[[k]][,i,j]) , min(post[[k]][,i,j]) )
                        result[r,4] <- the_spark
                        result[r,5] <- max( post[[k]][,i,j] )
                    }
                    rownames(result)[r] <- concat( names(post)[k] , "[" , i , "," , j , "]" )
                    r <- r + 1
                }
            }
        }
    }
    if ( spark==FALSE )
        colnames(result)[3:4] <- c(paste("lower", prob), paste("upper", prob))
    result
}

setGeneric("precis",
function( object , depth=1 , pars , prob=0.89 , digits=2 , sort=NULL , decreasing=FALSE , ... ) {   
    if ( "CmdStanFit" %in% class(object) ) {
        if ( missing(pars) ) pars <- NULL
        precis_cmdstanfit( object , depth=depth , pars=pars , prob=prob , digits=digits , sort=sort, decreasing=decreasing , ... )
    } else
        new( "precis" , as.data.frame(object) , digits=digits )
}
)

# function that handles depth and sorting
precis_format <- function( result , depth , sort , decreasing ) {
    # deal with depth
    if ( depth==1 ) {
        hits <- regexpr("]",rownames(result),fixed=TRUE)
        hits_idx <- which( hits > -1 )
        if ( length(hits_idx)>0 ) {
            result <- result[-hits_idx,]
            message( paste( length(hits_idx) , "vector or matrix parameters hidden. Use depth=2 to show them." ) )
        }
    }
    if ( depth==2 ) {
        hits <- regexpr(",",rownames(result),fixed=TRUE)
        hits_idx <- which( hits > -1 )
        if ( length(hits_idx)>0 ) {
            result <- result[-hits_idx,]
            message( paste( length(hits_idx) , "matrix parameters hidden. Use depth=3 to show them." ) )
        }
    }

    # sort
    if ( !is.null(sort) ) {
        o <- order( result[,sort] , decreasing=decreasing )
        result <- result[o,]
    }

    # label Rhat with version
    rhat_col <- which( colnames(result)=="Rhat" )
    if ( !is.null(rhat_col) ) {
        colnames(result)[rhat_col] <- "Rhat4"
    }

    return(result)
}

setMethod("precis", "numeric", 
function( object , depth=1 , pars , prob=0.89 , digits=2 , sort=NULL , decreasing=FALSE , ... ) {
    oname <- deparse( match.call()[[2]] )
    df <- list()
    df[[oname]] <- object
    precis( df , prob=prob , ... )
} )

setMethod("precis", "data.frame", 
function( object , depth=1 , pars , prob=0.89 , digits=2 , sort=NULL , decreasing=FALSE , hist=TRUE , ... ) {
    plo <- (1-prob)/2
    phi <- 1 - plo
    # replace any character or factor columns with NA numeric columns
    # histospark will detect all NA and return blank line
    for ( i in 1:ncol(object) ) {
        if ( class(object[[i]]) %in% c("factor","character") ) {
            object[[i]] <- as.numeric( rep( NA , nrow(object) ) )
        }
    }
    if ( hist==TRUE ) {
        result <- data.frame(
            mean = apply(object,2,mean,na.rm=TRUE),
            sd = apply(object,2,sd,na.rm=TRUE),
            lo = apply(object,2,quantile,na.rm=TRUE,probs=plo),
            hi = apply(object,2,quantile,na.rm=TRUE,probs=phi),
            histogram = apply(object,2,histospark),
            stringsAsFactors = FALSE
        )
    } else {
        # no unicode histogram
        result <- data.frame(
            mean = apply(object,2,mean,na.rm=TRUE),
            sd = apply(object,2,sd,na.rm=TRUE),
            lo = apply(object,2,quantile,na.rm=TRUE,probs=plo),
            hi = apply(object,2,quantile,na.rm=TRUE,probs=phi),
            stringsAsFactors = FALSE
        )
    }

    colnames(result)[3:4] <- paste( c( plo , phi )*100 , "%" , sep="" )

    result <- precis_format( result , depth , sort , decreasing )

    has_source <- !is.null(attr(object,"source"))
    header_string <- concat( "'data.frame': " , nrow(object) , " obs. of ", ncol(object) , " variables:" )
    if ( has_source )
        header_string <- attr(object,"source")
    attr(result,"header") <- header_string

    return( new( "precis" , result , digits=digits ) )
})

setMethod("precis", "list", 
function( object , depth=1 , pars , prob=0.89 , digits=2 , sort=NULL , decreasing=FALSE , hist=TRUE , ... ) {
    # coerce to data frame and format row names for vectors/matrices to [] style
    result <- as.data.frame( object , stringsAsFactors = FALSE )
    # since data frame conversion vectorizes matrices, need to treat each variable
    for ( i in 1:length(object) ) {
        # check dimension and process names when > 1
        n <- length(dim(object[[i]]))
        if ( n > 1 ) {
            dims <- dim(object[[i]])
            idx <- grep( concat("^",names(object)[i],".") , names(result) )
            if ( n==2 ) {
                # vector
                new_names <- paste( names(object)[i] , "[" , 1:dims[2] , "]" , sep="" )
                names(result)[idx] <- new_names
            }#2
            if ( n==3 ) {
                # matrix
                new_names <- paste( names(object)[i] , "[" , rep(1:dims[2],each=dims[3]) , "," , rep(1:dims[3],times=dims[2]) , "]" , sep="" )
                names(result)[idx] <- new_names
            }#3
        }#n>1
    }#i
    # hand off to data frame method
    if ( !is.null(attr(object,"source")) )
        attr(result,"source") <- attr(object,"source")
    precis( result , depth , pars , prob , digits , sort, decreasing , hist=hist , ... )
})



# template function for processing objects with coef and vcov methods
xprecis_glm <- function( object , depth , pars , prob , digits , sort , decreasing , ... ) {
    plo <- (1-prob)/2
    phi <- 1 - plo
    z <- qnorm(phi)
    result <- data.frame(
        mean = coef(object),
        sd = se(object),
        lo = coef(object) - z*se(object),
        hi = coef(object) + z*se(object)
    )
    colnames(result)[3:4] <- paste( c( plo , phi )*100 , "%" , sep="" )

    result <- precis_format( result , depth , sort , decreasing )

    # obey pars list
    if ( !missing(pars) ) {
        # need to handle vector/matrix parameters
        # for each element in pars, add a copy with '[' on end
        # then use grep to white list any parameter that starts with element of pars
        pars_new <- paste( pars , "[" , sep="" )
        pars <- c( pars , pars_new )
        keep_idx <- rep(FALSE,nrow(result))
        for ( i in 1:length(pars) ) {
            keep_idx <- keep_idx | startsWith( rownames(result) , pars[i] )
        }
        result <- result[ keep_idx , ]
    }

    return( new( "precis" , result , digits=digits ) )
}

setMethod("precis", "map",
function( object , depth=1 , pars , prob=0.89 , digits=2 , sort=NULL , decreasing=FALSE , ... ) {
    # just pass to glm method
    xprecis_glm( object , depth , pars , prob , digits , sort , decreasing , ... )
})

setMethod("precis", "glm",
function( object , depth=1 , pars , prob=0.89 , digits=2 , sort=NULL , decreasing=FALSE , ... ) {
    # just pass to glm method
    xprecis_glm( object , depth , pars , prob , digits , sort , decreasing , ... )
})

setMethod("precis", "lm",
function( object , depth=1 , pars , prob=0.89 , digits=2 , sort=NULL , decreasing=FALSE , ... ) {
    # just pass to glm method
    xprecis_glm( object , depth , pars , prob , digits , sort , decreasing , ... )
})

setMethod("precis", "map2stan",
function( object , depth=1 , pars , prob=0.89 , digits=2 , sort=NULL , decreasing=FALSE , ... ) {
    low <- (1-prob)/2
    upp <- 1-low
    result <- summary(object@stanfit,pars=pars,probs=c(low,upp))$summary[,c(1,3:7)]
    result <- as.data.frame( result )

    idx <- which( rownames(result) %in% c("dev","lp__") )
    idx2 <- grep( "log_lik[" , rownames(result) , fixed=TRUE )
    if ( length(idx2)>0 ) idx <- c( idx , idx2 )
    if ( length(idx)>0 ) {
        # remove dev and lp__ and log_lik from table
        result <- result[ -idx , ]
    }

    result <- precis_format( result , depth , sort , decreasing )

    return( new( "precis" , result , digits=digits ) )
})

setMethod("precis", "stanfit",
function( object , depth=1 , pars , prob=0.89 , digits=2 , sort=NULL , decreasing=FALSE , ... ) {
    # assume stan for moment
    low <- (1-prob)/2
    upp <- 1-low
    result <- summary(object,pars=pars,probs=c(low,upp))$summary[,c(1,3:7)]
    result <- as.data.frame( result )

    idx <- which( rownames(result) %in% c("dev","lp__") )
    idx2 <- grep( "log_lik[" , rownames(result) , fixed=TRUE )
    if ( length(idx2)>0 ) idx <- c( idx , idx2 )
    if ( length(idx)>0 ) {
        # remove dev and lp__ and log_lik from table
        result <- result[ -idx , ]
    }

    result <- precis_format( result , depth , sort , decreasing )

    return( new( "precis" , result , digits=digits ) )
})
# precis2(mfit,pars=parlist)

# cmdstanr doesn't export CmdStanFit class, so need to do this through default method
#setMethod("precis", "CmdStanFit",
precis_cmdstanfit <- function( object , depth=1 , pars , prob=0.89 , digits=2 , sort=NULL , decreasing=FALSE , ... ) {
    low <- (1-prob)/2
    upp <- 1-low
    # result <- summary(object,variables=pars,probs=c(low,upp))$summary[,c(1,3:7)]
    result <- object$summary( variables=pars , "mean" , "sd" , ~quantile(.x, probs = c(low, upp)) , "rhat" , "ess_bulk" )
    result <- as.data.frame( result )

    # take variable column and convert to row names
    rownames(result) <- result$variable
    result$variable <- NULL

    idx <- which( rownames(result) %in% c("dev","lp__") )
    idx2 <- grep( "log_lik[" , rownames(result) , fixed=TRUE )
    if ( length(idx2)>0 ) idx <- c( idx , idx2 )
    if ( length(idx)>0 ) {
        # remove dev and lp__ and log_lik from table
        result <- result[ -idx , ]
    }

    result <- precis_format( result , depth , sort , decreasing )

    return( new( "precis" , result , digits=digits ) )
}

# old version 1 method
precisx <- function( model , depth=1 , pars , ci=TRUE , prob=0.89 , corr=FALSE , digits=2 , warn=TRUE , spark=FALSE ) {
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
        message( paste("No handler found for model of class",the.class) )
        return(invisible())
    }
    # format
    fname <- "Mean"
    # capitalize first letter
    fname <- concat( toupper(substring(fname,1,1)) , substring(fname,2) )
    result <- data.frame( est=est , se=se )
    colnames(result) <- c( fname ,"StdDev")
    if ( ci==TRUE ) {
        ci <- confint_quad( est=est , se=se , prob=prob )
        if ( the.class=="data.frame" ) {
            # HPDI from samples
            ci <- t( apply( model , 2 , HPDI , prob=prob ) )
        }
        result <- cbind( result , ci )
        if ( the.class=="map2stan" ) {
            # HPDI from samples
            post <- extract.samples(model)
            result <- postlistprecis( post , prob=prob , spark=spark )
        }
        if ( the.class=="stanfit" ) {
            # HPDI from samples
            post <- extract.samples(model)
            post[['lp__']] <- NULL
            result <- postlistprecis( post , prob=prob , spark=spark )
        }
    }
    if ( the.class=="map2stan" | the.class=="stanfit" ) {
        # add n_eff to result
        #require(rstan)
        if ( the.class=="map2stan" )
            the_summary <- summary( model@stanfit )$summary
        else
            the_summary <- summary( model )$summary
        n_eff <- the_summary[,'n_eff']
        n_eff <- n_eff[ -which(names(n_eff)=="lp__") ]
        Rhat <- the_summary[,'Rhat']
        Rhat <- Rhat[ -which(names(Rhat)=="lp__") ]
        if ( the.class=="map2stan" ) {
            n_eff <- n_eff[ -which(names(n_eff)=="dev") ]
            Rhat <- Rhat[ -which(names(Rhat)=="dev") ]
        }
        result <- cbind( result , n_eff , Rhat )
        
        # check divergent iterations
        nd <- divergent(model)
        if ( nd > 0 & warn==TRUE ) {
            warning( concat("There were ",nd," divergent iterations during sampling.\nCheck the chains (trace plots, n_eff, Rhat) carefully to ensure they are valid.") )
        }
    }
    if ( corr==TRUE ) {
        result <- cbind( result , Rho )
    }
    
    #if ( type.s==TRUE )
    #    result[,"Pr(S)"] <- format.pval( type.s( est , se ) )
        
    if ( precis.whitelist$vcov.method[ precis.whitelist$class==the.class ]=="vcov.VarCorr" ) {
        message( "Quadratic approximation (standard errors) unreliable for variance components. Use MCMC to estimate precision of variance components." )
    }
    
    # deal with depth
    if ( depth==1 ) {
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
rrformat <- function( matrix , digits=2 , width=7 ) {
    if ( length(digits)==1 ) digits <- rep(digits,nrow(matrix))
    result <- matrix
    for ( i in 1:nrow(matrix) ) {
        result[i,] <- format( round(matrix[i,],digits[i]) , width=width )
    }
    result
}

