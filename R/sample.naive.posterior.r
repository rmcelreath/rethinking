# sample naive posterior

# convenience wrapper for sampling from naive posterior into a data frame
# also capable of averaging over models, if given a list() as first argument
sample.naive.posterior <- function( model , n=10000 , clean.names=TRUE , method="AICc" , nobs=0 , add.names=FALSE , fill.na=0 , verbose=FALSE ) {
    require(MASS)
    require(bbmle)
    # need own BIC, as one in stats doesn't allow nobs
    getdf <- function(x) {
        if (!is.null(df <- attr(x, "df"))) 
            return(df)
        else if (!is.null(df <- attr(logLik(x), "df"))) 
            return(df)
    }
    myBIC <- function(x,nobs) {
        k <- getdf(x)
        as.numeric( -2*logLik(x) + log(nobs)*k )
    }
    if ( class(model)[1]=="list" ) {
        # list of models passed
        valid.methods <- c("AIC","AICc","BIC")
        if ( !any(method==valid.methods) ) {
            cat( paste("Unknown model averaging method:",method) )
            return()
        }
        if ( length( model ) < 2 ) {
            return( sample.naive.posterior( model[[1]] , n=n , method=method , nobs=nobs ) )
        }
        if ( method=="AIC" ) factors <- sapply( model , AIC )
        if ( nobs==0 ) {
            if ( method=="AICc" ) factors <- sapply( model , AICc )
            if ( method=="BIC" ) factors <- sapply( model , BIC )
        } else {
            if ( method=="AICc" ) factors <- sapply( model , function(z) AICc(z,nobs=nobs) )
            if ( method=="BIC" ) factors <- sapply( model , function(z) myBIC(z,nobs=nobs) )
        }
        factors <- factors - min( factors )
        post <- exp( -1/2 * factors )/sum( exp( -1/2 * factors ) )
        sim.post <- vector( "list" , length(model) )
        f.zeros <- FALSE
        nn <- round(n*post)
        if ( verbose ) print(nn)
        for ( i in 1:length(model) ) {
            if ( nn[i]==0 ) {
                f.zeros <- TRUE
                sim.post[[i]] <- sample.naive.posterior( model[[i]] , n=2 ) 
                # 2 appears to be minimum to get columns to populate from mvrnorm
                # need column names from these zero samples models
            } else {
                sim.post[[i]] <- sample.naive.posterior( model[[i]] , n=max(nn[i],2) )
            }
        }
        if ( f.zeros==TRUE ) {
            warning( "One or more models produced zero samples, because of very low posterior probability." )
        }
        # new algorithm
        if ( TRUE ) {
            # build list of unique parameter names
            par.names <- sapply( sim.post , function(i) colnames( i ) )
            upar.names <- unique( unlist( par.names ) )
            # build averaged posterior
            post.avg <- matrix( NA , nrow=sum(nn) , ncol=length(upar.names) )
            colnames(post.avg) <- upar.names
            # insert samples
            current.row <- 1
            for ( i in 1:length(model) ) {
                if ( nn[i] > 0 ) {
                    start.row <- current.row
                    end.row <- current.row + nn[i] - 1
                    for ( j in colnames(sim.post[[i]]) ) {
                        post.avg[ start.row:end.row , j ] <- sim.post[[i]][ 1:nn[i] , j ]
                    }
                    current.row <- end.row + 1
                }
            }
            post.avg <- as.data.frame( post.avg )
        } else {
        # old algorithm
            post.avg <- sim.post[[1]]
            # post.avg <- merge( sim.post[[1]] , sim.post[[2]] , all=TRUE , sort=FALSE )
            if ( length( model ) > 1 ) {
                for ( i in 2:length(model) ) {
                    # if ( nn[i] > 0 ) 
                    post.avg <- merge( post.avg , sim.post[[i]] , all=TRUE , sort=FALSE )
                    if ( nn[i] == 0 ) {
                        nr <- nrow(post.avg)
                        rcut <- (nr-1):nr
                        post.avg <- post.avg[ -rcut , ]
                    }
                }
                if ( nn[1] == 0 ) {
                    post.avg <- post.avg[ -(1:2) , ]
                }
            } # end old algorithm
        }
        if ( !is.logical(fill.na) )
            post.avg[is.na(post.avg)] <- fill.na
        if ( add.names==TRUE ) {
            mnames <- match.call()
            mnames <- as.character(mnames[[2]])[2:(length(model)+1)]
            rnames <- rep( mnames , times=nn )
            post.avg$model <- rnames
        }
        result <- post.avg
    } else {
        # single model passed
        mu <- 0
        if ( class(model)[1] %in% c("mer","bmer") ) {
            mu <- fixef(model)
        } else {
            mu <- xcoef(model)
        }
        result <- as.data.frame( mvrnorm( n=n , mu=mu , Sigma=vcov(model) ) )
        if ( clean.names==TRUE ) {
            # convert (Intercept) to Intercept
            for ( i in 1:ncol(result) ) {
                if ( colnames(result)[i] == "(Intercept)" ) {
                    colnames(result)[i] <- "Intercept"
                }
            }
        }
    }
    result
}
