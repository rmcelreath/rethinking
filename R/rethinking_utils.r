# various utility functions

# default quartz plot size for book: 3.5in by 4in, giving square plot for default margins
blank <- function(ex=1) quartz("myquartz",width=3.5*ex,height=4*ex)

# convenience function for choosing a csv file
choose.csv <- function( ... ) read.csv( file=file.choose() , ... )

# bound a list of real values on either side
fclip <- function( x , xmin=NULL , xmax=NULL ) {
    if ( !is.null(xmin) ) x <- ifelse( x < xmin , xmin , x )
    if ( !is.null(xmax) ) x <- ifelse( x > xmax , xmax , x )
    x
}

make.grid <- function( n ) {
    num.rows <- floor( sqrt(n) )
    num.cols <- ceiling(n/num.rows)
    c(num.rows,num.cols)
}

dens <- function( x , adj=0.1 , norm.comp=FALSE , main="" , show.HPDI=FALSE , show.zero=FALSE , rm.na=TRUE , ...) {
    the.class <- class(x)[1]
    if ( the.class=="data.frame" ) {
        # full posterior
        n <- ncol(x)
        cnames <- colnames(x)
        par( mfrow=make.grid(n) )
        for ( i in 1:n ) {
            dens( x[,i] , adj=adj , norm.comp=norm.comp , show.HPDI=show.HPDI , show.zero=TRUE , xlab=cnames[i] , ... )
        }
    } else {
        # vector
        if ( rm.na==TRUE ) x <- x[ !is.na(x) ]
        thed <- density(x,adjust=adj)
        if ( show.HPDI==FALSE )
            plot( thed , main=main , ... )
        else
            plot( thed , main=main , type="n" , ... )
        if ( show.HPDI != FALSE ) {
            hpd <- HPDI( x , prob=show.HPDI )
            i <- which( thed$x >= hpd[1] & thed$x <= hpd[2] )
            b.low <- min( i )
            b.hi <- max( i )
            for ( j in i ) {
                lines( c( thed$x[j] , thed$x[j] ) , c( 0 , thed$y[j] ) , col=col.alpha("slateblue",0.2) )
            }
            lines( thed$x , thed$y )
        }
        if ( norm.comp==TRUE ) {
            mu <- mean(x)
            sigma <- sd(x)
            curve( dnorm( x , mu , sigma ) , col="white" , lwd=2 , add=TRUE )
            curve( dnorm( x , mu , sigma ) , add=TRUE )
        }
        if ( show.zero==TRUE ) {
            lines( c(0,0) , c(0,max(thed$y)*2) , lty=2 )
        }
    }
}

# just converts x,y,z lists of same length to a matrix for contour to plot
contour.xyz <- function( x , y , z , ... ) {
    ux <- unique(x)
    uy <- unique(y)
    n <- length(ux)
    m <- matrix( z , nrow=n , ncol=n )
    contour( ux , uy , m , ... )
}

# finds mode of a continuous density
chainmode <- function( chain , ... ) {
    dd <- density(chain , ...)
    dd$x[which.max(dd$y)]
}

# highest posterior density interval, sensu Box and Tiao
# requires coda library
HPDI <- function( samples , prob=0.95 ) {
    require(coda)
    class.samples <- class(samples)[1]
    coerce.list <- c( "numeric" , "matrix" , "data.frame" , "integer" )
    if ( class.samples %in% coerce.list ) {
        # single chain for single variable
        samples <- as.mcmc( samples )
    }
    x <- coda::HPDinterval( samples , prob=prob )
    result <- c( x[1] , x[2] )
    names(result) <- c(paste("lower",prob),paste("upper",prob))
    result
}

# percentile confidence interval
PCI <- function( samples , prob=0.95 ) {
    a <- (1-prob)/2
    quantile( samples , probs=c(a,1-a) )
}

# plot new y values on secondary y-axis on right side
plot2y <- function( ... , y2lab="2nd axis" , y2col=NULL ) {
    if ( is.null(y2col) ) y2col <- "black"
    par(new=TRUE)
    par(mar=c(5,4,4,5)+.1)
    plot( xaxt="n" , yaxt="n" , xlab="" , ylab="" , ... )
    axis( 4 , col=y2col , col.axis=y2col )
    mtext( y2lab , side=4 , line=3 , col=y2col )
}

# color utility functions
col.desat <- function( acol , amt=0.5 ) {
    acol <- col2rgb(acol)
    ahsv <- rgb2hsv(acol)
    ahsv[2] <- ahsv[2] * amt
    hsv( ahsv[1] , ahsv[2] , ahsv[3] )
}

col.alpha <- function( acol , alpha=0.2 ) {
    acol <- col2rgb(acol)
    acol <- rgb(acol[1]/255,acol[2]/255,acol[3]/255,alpha)
    acol
}

col.dist <- function( x , mu=0 , sd=1 , col="slateblue" ) {
    cols <- sapply( x , function(z) exp(-(z-mu)^2/sd) )
    cols <- sapply( cols , function(z) col.alpha(col,z) )
    cols
}

se <- function( model ) {
    sqrt( diag( vcov(model) ) )
}

# function for plotting "naive" posteriors and 95% confints
show.naive.posterior <- function( est , se , model=NULL , level=0.95 , xlab="estimate" , ylab="likelihood" , npts=1000 , ciy=NULL , show.density=TRUE , show.ci=TRUE , zero.lines=TRUE , label=NULL , cols=NULL , lwidths=NULL , ... ) {
    if ( !is.null(model) ) {
        f.found.class <- FALSE
        if ( class(model)=="lm" ) {
            f.found.class <- TRUE
            est <- coef(model)
            if ( is.null(label) ) label <- names(est)
            est <- as.vector(est)
            se <- as.vector(summary(model)$coefficients[,2])
        }
        if ( class(model)=="mle2" ) {
            f.found.class <- TRUE
            est <- coef(model)
            if ( is.null(label) ) label <- names(est)
            est <- as.vector(est)
            se <- as.vector(summary(model)@coef[,2])
        }
        if ( class(model)[1]=="mer" ) {
            f.found.class <- TRUE
            est <- fixef(model)
            if ( is.null(label) ) label <- names(est)
            est <- as.vector(est)
            se <- as.vector( sqrt( diag( vcov(model) ) ) )
        }
        if ( f.found.class==FALSE ) {
            return( paste("Could not find handler for model of class",class(model)) )
        }
    }
    minx <- min( est - se*3 )
    maxx <- max( est + se*3 )
    y <- matrix( 0 , nrow=length(est) , ncol=npts )
    ci <- matrix( 0 , nrow=length(est) , ncol=2 )
    tails <- (1-level)/2
    tails <- c( tails , 1-tails )
    x <- seq( from=minx , to=maxx , length=npts )
    for ( i in 1:nrow(y) ) {
        y[i,] <- dnorm( x , est[i] , se[i] )
        ci[i,] <- qnorm( tails , est[i] , se[i] )
    }
    maxy <- max( y )
    if ( is.null(ciy) ) ciy <- -maxy/10
    if ( show.ci==FALSE ) ciy <- 0
    yaxistype <-"s"
    if ( show.density==FALSE ) {
        maxy <- 0
        yaxistype <- "n"
        ylab=""
    }
    plot( 0 , 0 , type="n" , xlab=xlab , ylab=ylab , ylim=c(ciy,maxy) , xlim=c(minx,maxx) , yaxp=c(0,maxy,4) , yaxt=yaxistype , ... )
    if ( zero.lines==TRUE ) {
        if ( show.density==TRUE ) 
            lines( c(minx-abs(minx),maxx+abs(maxx)) , c(0,0) , lty=3 )
        lines( c(0,0) , c(-1,maxy*2) , lty=3 )
    }
    if ( is.null(cols) ) cols <- rep( "black" , length(est) )
    if ( is.null(lwidths) ) lwidths <- rep( 1 , length(est) )
    if ( show.density==TRUE ) {
        for ( i in 1:nrow(y) ) {
            lines( x , y[i,] , col=cols[i] , lwd=lwidths[i] , ... )
        }
    }
    yoff.ci <- rep(0,nrow(ci))
    if ( show.ci==TRUE ) {
        for ( i in 1:nrow(ci) ) {
            yoff.ci[i] <- ciy/(nrow(ci)+1)*i
            lines( ci[i,] , rep( yoff.ci[i] , 2 ) , col=cols[i] , lwd=lwidths[i] , ... )
            points( est[i] , yoff.ci[i] , col=cols[i] , lwd=lwidths[i] , ... )
        }
    }
    if ( !is.null(label) ) {
        xloc <- est
        yloc <- dnorm( xloc , xloc , se )
        yoff <- ifelse( yloc==max(yloc) , ciy/2 , -ciy/2 )
        yloc <- yloc + yoff
        if ( show.density==FALSE ) {
            yloc <- yoff.ci + 0.01
        }
        for ( i in 1:length(label) ) {
            text( xloc[i] , yloc[i] , labels=label[i] , col=cols[i] )
        }
    }
}

# simple histogram
simplehist <- function( x , ylab="Frequency" , xlab="Count" , ycounts=TRUE , adjust=1 , lcol="black" , bins=NULL , show.counts=0 , xlim=NULL , ylim=NULL , ... ) {
    # first, check if integers only or continuous.
    freqs <- {}
    x2 <- round(x)
    iflag <- FALSE
    if ( all(x2==x) | !is.null(bins) ) {
        # integer dist, so plot each value
        if ( is.null(bins) ) {
            bins <- min(0,x):max(x)
        }
        freqs <- sapply( bins , function(z) length(x[x==z])/sum(x) )
        if ( ycounts==TRUE ) freqs <- sapply( bins , function(z) length( x[ as.character(x)==as.character(z) ] ) )
        iflag <- TRUE
    } else {
        # continuous dist (fractional values detected)
    }
    # plot frame
    if ( iflag==TRUE ) {
        # integers
        if ( is.null(ylim) )
            ylim <- c( 0 , max(freqs) )
        if ( is.null(xlim) )
            xlim <- c(min(bins),max(bins))
        plot( 0 ~ 0 , col="white" , xlim=xlim , ylim=ylim , ylab=ylab , xlab=xlab , ... )
        for ( i in 1:length(bins) ) {
            if ( freqs[i] > 0 )
                lines( c(bins[i],bins[i]) , c(0,freqs[i]) , col=lcol , ... )
        }
        # finally, show text counts for all bars with counts < show.counts
        if ( show.counts > 0 ) {
            for ( i in 1:length(bins) ) {
                if ( freqs[i] < show.counts )
                    text( bins[i] , freqs[i] , labels=freqs[i] , pos=3 )
            }
        }
    } else {
        # continuous
       if ( ylab=="Frequency" ) ylab <- "Density"
       if ( xlab=="Count" ) xlab <- "Value"
       plot( density( x , adjust=adjust ) , ylab=ylab , xlab=xlab , ... )
    }
    
}

# quadratic estimate confidence intervals from means and standard errors
confint.quad <- function( model=NULL , est , se , level=0.95 ) {
    if ( !is.null(model) ) {
        found.class <- FALSE
        if ( class(model)=="lm" ) {
            est <- coef(model)
            se <- summary(model)$coef[,2]
            found.class <- TRUE
        }
        if ( class(model)=="mle2" ) {
            est <- coef(model)
            se <- summary(model)@coef[,2]
            found.class <- TRUE
        }
        if ( found.class==FALSE ) {
            return( paste("Cannot find handler for model of class",class(model)) )
        }
    }
    n <- length(est)
    mat <- matrix(c(rep(-1,n),rep(1,n)),nrow=n)
    p <- (1-level)/2
    z <- -qnorm( p )
    ci <- est + mat * ( se * z )
    rownames(ci) <- names(est)
    lowlab <- paste( format( p*100 , nsmall=1 ) , "%" , sep="" )
    hilab <- paste( format( (1-p)*100 , nsmall=1 ) , "%" , sep="" )
    colnames(ci) <- c( lowlab , hilab )
    ci
}

# Type S error function
type.s <- function( est , se , posterior=NULL ) {
    n <- length(est)
    if ( is.null(posterior) ) {
        pr.s <- pnorm( rep(0,n) , mean=est , sd=se )
        pr.s <- ifelse( pr.s > 0.5 , 1-pr.s , pr.s )
    } else {
        # compute from posterior rather than standard error
    }
    pr.s
}

# my model summary function, précis
precis.whitelist <- data.frame( 
    class=c("lm","glm","mle2","mer","bmer","polr","data.frame","clmm","clmm2","list","stanfit","lmerMod","glmerMod") , 
    coef.method=c("coef","coef","coef","fixef.plus","fixef.plus","polr","chain","coef","coef","mcarray","stanfit","fixef.plus","fixef.plus") , 
    vcov.method=c("vcov","vcov","vcov","vcov.VarCorr","vcov.VarCorr","vcov","chain","vcov","vcov","mcarray","stanfit","vcov.VarCorr","vcov.VarCorr") ,
    nobs.method=c("nobs","nobs","mle2","mer","mer","nobs","chain","nobs","nobs","chain","stanfit","mer","mer")
)

# precis class definition and show method
setClass( "precis" , representation( output="data.frame" , digits="numeric" ) )
precis.show <- function( object ) {
    print( round( object@output , object@digits ) )
}
setMethod( "show" , signature=(x="precis") , function(object) precis.show(object) )

precis <- function( model , type.s=FALSE , ci=TRUE , level=0.95 , digits=2 , warn=TRUE ) {
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
    if ( type.s==TRUE )
        result[,"Pr(S)"] <- format.pval( type.s( est , se ) )
    if ( precis.whitelist$vcov.method[ precis.whitelist$class==the.class ]=="vcov.VarCorr" ) {
        message( "Quadratic approximation (standard errors) unreliable for variance components. Use MCMC to estimate precision of variance components." )
    }
    new( "precis" , output=result , digits=digits )
}

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

mcmcpairs <- function( posterior , cex=0.3 , pch=16 , col=col.alpha("slateblue",0.2) , n=1000 , ... ) {
    panel.dens <- function(x, ...) {
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(usr[1:2], 0, 1.5) )
        h <- density(x,adj=0.5)
        y <- h$y
        y <- y/max(y)
        abline( v=0 , col="gray" , lwd=0.5 )
        lines( h$x , y )
    }
    panel.2d <- function( x , y , ... ) {
        i <- sample( 1:length(x) , size=n )
        abline( v=0 , col="gray" , lwd=0.5 )
        abline( h=0 , col="gray" , lwd=0.5 )
        points( x[i] , y[i] , ... )
    }
    panel.cor <- function( x , y , ... ) {
        k <- cor( x , y )
        cx <- sum(range(x))/2
        cy <- sum(range(y))/2
        text( cx , cy , round(k,2) , cex=2*exp(abs(k))/exp(1) )
    }
    pairs( posterior , cex=cex , pch=pch , col=col , upper.panel=panel.2d , lower.panel=panel.cor , diag.panel=panel.dens , ... )
}

# convert to sigma from tau
tau2sig <- function( tau ) 1/abs(tau)

# Gaussian density that takes tau instead of sigma
dnormtau <- function( x , mean=0 , tau=1 , log=FALSE ) dnorm( x , mean=mean , sd=tau2sig( tau ) , log=log )

rnormtau <- function( n , mean=0 , tau=1 ) rnorm( n=n , mean=mean , sd=tau2sig(tau) )

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

# coeftab class definition and show method
setClass( "coeftab" , representation( coefs="matrix" , nobs="numeric" , AIC="numeric" , digits="numeric" , width="numeric" ) )
coeftab.show <- function( object ) {
    result <- object@coefs
    if ( !is.null(object@nobs) ) {
        result <- rbind( result , object@nobs )
        rownames(result)[ nrow(result) ] <- "nobs"
    }
    coefs <- rrformat( result , digits=object@digits , width=object@width )
    print( coefs , quote=FALSE , justify="right" )
}
setMethod( "show" , signature=(x="coeftab") , function(object) coeftab.show(object) )

coeftab <- function( ... , se=FALSE , se.inside=FALSE , nobs=TRUE , digits=2 , width=7 , rotate=FALSE , compare=FALSE ) {
    
    # se=TRUE outputs standard errors
    # se.inside=TRUE prints standard errors in parentheses in same column as estimates
    
    if ( se.inside==TRUE ) se <- TRUE
    
    # retrieve list of models
    L <- list(...)
    if ( is.list(L[[1]]) && length(L)==1 )
        L <- L[[1]]
    
    # retrieve model names from function call
    mnames <- match.call()
    mnames <- as.character(mnames)[2:(length(L)+1)]
    
    # count number of unique parameters
    param.names <- {}
    for ( i in 1:length(L) ) {
        c.names <- names( xcoef( L[[i]] ) )
        param.names <- unique( c( param.names , c.names ) )
    }
    # columns for standard errors
    if ( se==TRUE && se.inside==FALSE ) {
        for ( i in 1:length(L) ) {
            kse.names <- paste( names( xcoef( L[[i]] ) ) , ".se" , sep="" )
            param.names <- unique( c( param.names , kse.names ) )
        }
    }
    
    # make empty table
    nk <- length(param.names)
    d <- matrix( NA , ncol=nk )
    d <- data.frame(d)
    colnames(d) <- c( param.names )
    
    # loop over models and insert values
    for ( i in 1:length(L) ) {
        klist <- xcoef( L[[i]] )
        for ( j in 1:length(klist) ) {
            d[i,][ names( klist[j] ) ] <- as.numeric( round( klist[j] , digits ) )
        }
    }
    # insert standard errors
    if ( se==TRUE ) {
        for ( i in 1:length(L) ) {
            kse <- xse( L[[i]] )
            names(kse) <- names( xcoef( L[[i]] ) )
            for ( j in 1:length(kse) ) {
                if ( se.inside==FALSE )
                    # own column
                    d[i,][ paste(names(kse)[j],".se",sep="") ] <- as.numeric( round(kse[j],digits) )
                else
                    # combine with estimate
                    d[i,][ names(kse)[j] ] <- paste( formatC( as.real(d[i,][ names(kse)[j] ]) , digits=digits ) , " (" , formatC( as.real( kse[j] ) , digits=digits ) , ")" , sep="" )
            }
        }
    }
    
    # add model names to rows
    rownames(d) <- mnames
    
    # formatting for parenthetical standard errors
    if ( se.inside==TRUE && se==TRUE ) {
        comment(d) <- "Values in parentheses are quadratic estimate standard errors."
        colnames(d) <- paste( colnames(d) , "(se)" )
        for ( i in 1:nrow(d) ) {
            for ( j in 1:ncol(d) ) {
                d[i,j] <- ifelse( is.na(d[i,j]) , "" , d[i,j] )
            }
        }
    }
    
    # add nobs
    if ( nobs==TRUE ) {
        nobs <- sapply( L , xnobs )
    } else {
        nobs <- 0
    }
    
    # add AICc and weights
    if ( compare==TRUE ) {
        require(bbmle)
        d$AICc <- sapply( 1:length(L) , function(z) AICc( L[[z]] , nobs=d$nobs[z] ) )
        deltAICc <- d$AICc - min(d$AICc) 
        wAIC <- exp( -0.5*deltAICc )
        d$weight <- wAIC / sum( wAIC )
    }
    
    # return table
    if ( rotate==FALSE ) d <- t(d) # models along top is default
    new( "coeftab" , coefs=d , nobs=nobs , digits=digits , width=width )
}


######
# take a fit lm object and return a fit mle2 object with same model, data
# uses estimates from lm as starting values for mle2
# uses do.call()
lm.to.mle2 <- function( model , data=NULL , tau=TRUE , skip.hessian=FALSE ) {
    if ( class(model)[1] == "formula" ) {
        # passed formula directly, so need to fit it ourselves
        model <- lm( model , data=data )
    }
    if ( class(model)[1] != "lm" ) return( "Not a model of class lm." )
    if ( is.null(data) ) return( "Must specify data frame." )
    require(bbmle)
    # extract data frame
    dm <- model$model
    # get factors matrix
    fm <- attr( attr( model$model , "terms" ) , "factors" )
    # get names of variables
    xnames <- colnames(fm)
    yname <- rownames(fm)[1]
    # build vector of additive terms for formula
    k <- length(xnames)
    terms <- rep("",k)
    coefnames <- paste( "b" , 1:k , sep="" )
    for ( i in 1:k ) {
        aterm <- ifelse( fm[,i]==1 , colnames(dm) , "" )
        aterm <- aterm[ aterm != "" ]
        acoef <- coefnames[i]
        aterm <- paste( c(acoef,aterm) , collapse="*" )
        terms[i] <- aterm
    }
    mu <- paste( c("a",terms) , collapse=" + " )
    # make start list
    startvals <- coef(model)
    sigmastart <- sqrt( var(resid(model))*(nrow(dm)-1)/(nrow(dm)-k-1) )
    sigmaname <- "sigma"
    sigmaform <- "sigma"
    if ( tau==TRUE ) {
        sigmaname <- "tau"
        sigmastart <- 1/sigmastart
        sigmaform <- "1/abs(tau)"
    }
    thestart <- paste( c("a",coefnames,sigmaname) , c( startvals , sigmastart ) , sep="=" )
    dostart <- paste( thestart , collapse=" , " )
    dostart <- paste( "list(" , dostart , ")" , sep="" )
    dostart <- eval( parse( text=dostart ) )
    # build model formula
    theform <- paste( yname , "~ dnorm( mean=" , mu , " , sd=" , sigmaform , ")" )
    theformula <- formula( theform )
    # put it all together
    holdcall <- match.call()
    fit.m <- do.call( mle2 , list(minuslogl=theformula , data=data , start=dostart,skip.hessian=skip.hessian) )
    fit.m@call <- holdcall
    fit.m@call.orig <- holdcall
    return( fit.m )
}

# m2 <- lm.to.mle2( lm( dist ~ speed , data=d ) , data=d )

# shannon information
#H <- function( p ) {
#    u <- ifelse( p==0 , 0 , log(p) )
#    -sum( p * u )
#}
# conditional entropy - uncertainty in X knowing Y
#HXY <- function( p ) {
#    pY <- apply( p , 2 , sum )
#    h <- 0
#    for ( i in 1:nrow(p) ) {
#        for ( j in 1:ncol(p) ) {
#            h <- h + p[i,j] * log( pY[j] / p[i,j] )
#        }
#    }
#    h
#}
# K-L divergence
# p is true function
# q is approximating function
#DKL <- function( p , q ) sum( p * log( p / q ) )

# ordered categorical density functions
logistic <- function( x ) {
    p <- exp( x ) / ( 1 + exp( x ) )
    p <- ifelse( x==Inf , 1 , p )
    p
}

pordlogit <- function( x , a , phi , log=FALSE ) {
    a <- c( as.numeric(a) , Inf )
    if ( length(phi) == 1 ) {
        p <- logistic( a[x] - phi )
    } else {
        p <- matrix( NA , ncol=length(x) , nrow=length(phi) )
        for ( i in 1:length(phi) ) {
            p[i,] <- logistic( a[x] - phi[i] )
        }
    }
    if ( log==TRUE ) p <- log(p)
    p
}

dordlogit <- function( x , a , phi , log=FALSE ) {
    a <- c( as.numeric(a) , Inf )
    p <- logistic( a[x] - phi )
    na <- c( -Inf , a )
    np <- logistic( na[x] - phi )
    p <- p - np
    if ( log==TRUE ) p <- log(p)
    p
}

rordlogit <- function( n , a , phi=0 ) {
    a <- c( as.numeric(a) , Inf )
    k <- 1:length(a)
    p <- dordlogit( k , a=a , phi=phi , log=FALSE )
    y <- sample( k , size=n , replace=TRUE , prob=p )
    y
}

# compare class definition and show method
setClass( "compareIC" , representation( output="data.frame" ) )
compare.show <- function( object ) {
    print( round( object@output , 2 ) )
}
setMethod( "show" , signature=(x="compareIC") , function(object) compare.show(object) )

# AICc/BIC model comparison table
compare <- function( ... , nobs=NULL , sort="AICc" , delta=TRUE , digits=4 ) {
    require(bbmle)
    
    if ( is.null(nobs) ) {
        stop( "Must specify number of observations (nobs)." )
    }
    
    getdf <- function(x) {
        if (!is.null(df <- attr(x, "df"))) 
            return(df)
        else if (!is.null(df <- attr(logLik(x), "df"))) 
            return(df)
    }
    
    # need own BIC, as one in stats doesn't allow nobs
    myBIC <- function(x,nobs) {
        k <- getdf(x)
        as.numeric( -2*logLik(x) + log(nobs)*k )
    }
    
    # retrieve list of models
    L <- list(...)
    if ( is.list(L[[1]]) && length(L)==1 )
        L <- L[[1]]
    
    # retrieve model names from function call
    mnames <- match.call()
    mnames <- as.character(mnames)[2:(length(L)+1)]
    
    AICc.list <- sapply( L , function(z) AICc( z , nobs=nobs ) )
    BIC.list <- sapply( L , function(z) myBIC( z , nobs=nobs ) )
    dAICc <- AICc.list - min( AICc.list )
    dBIC <- BIC.list - min( BIC.list )
    
    post.AICc <- exp( -0.5*dAICc ) / sum( exp(-0.5*dAICc) )
    post.BIC <- exp( -0.5*dBIC ) / sum( exp(-0.5*dBIC) )
    
    #post.AICc <- format.pval( post.AICc , digits=digits )
    #post.BIC <- format.pval( post.BIC , digits=digits )
    
    k <- sapply( L , getdf )
    
    result <- data.frame( k=k , AICc=AICc.list , BIC=BIC.list , w.AICc=post.AICc , w.BIC=post.BIC )
    if ( delta==TRUE ) {
        r2 <- data.frame( dAICc=dAICc , dBIC=dBIC )
        result <- cbind( result , r2 )
    }
    rownames( result ) <- mnames
    
    if ( !is.null(sort) ) {
        result <- result[ order( result$AICc ) , ]
    }
    
    new( "compareIC" , output=result )
}

# convenience interface to glmer (lme4) that always uses REML=FALSE
glmm <- function( ... , family , REML=FALSE ) {
    require(lme4)
    if ( missing(family) ) {
        result <- lmer( ... , REML=REML )
    } else {
        result <- glmer( ... , family=family )
    }
    result
}

bglmm <- function( ... , REML=FALSE ) {
    require(lme4)
    require(blme)
    blmer( ... , REML=REML )
}

# timing functions

progbar <- function( current , min=0 , max=100 , starttime , update.interval=100 , show.rate=FALSE ) {
    z <- current/update.interval
    if ( floor(z)!=z ) return()
    progress <- current - min
    time.elapsed <- Sys.time() - starttime
    elapsed.units <- attr(time.elapsed,"units")
    f <- 1
    if ( elapsed.units == "mins" ) f <- 60
    if ( elapsed.units == "hours" ) f <- 60*60
    if ( elapsed.units == "days" ) f <- 60*60*24
    rate <- progress / ( as.numeric(time.elapsed) * f )
    eta <- (max - current)/rate
    etamins <- floor(eta/60)
    etahours <- floor(etamins/60)
    etasecs <- floor(eta - etamins*60)
    if ( etahours > 0 ) etamins <- etamins - etahours*60
    percentdone <- floor(current/max * 100)
    ETAstring <- paste( etahours , "h " , sep="" )
    if ( etahours==0 & etamins>0 ) ETAstring <- paste( etamins , "m " , sep="" )
    if ( etahours==0 & etamins==0 ) ETAstring <- paste( etasecs , "s " , sep="" )
    Ratestring <- ""
    if ( show.rate==TRUE ) Ratestring <- paste( " | Rate: " , round(rate,2) , " samples/sec" , sep="" )
    cat( "\r" )
    cat( paste( c( current , "/" , max , " (" , percentdone , "%) | Elapsed: " , round(time.elapsed,2) , elapsed.units , " | Remaining: " , ETAstring , Ratestring ) , collapse="" ) )
    #cat( paste( c( current , "/" , max , " (" , percentdone , "%) | Elapsed: " , round(time.elapsed,2) , elapsed.units , " | Remaining: " , etahours , "hrs " , etamins , "mins " , etasecs , "secs | Rate: " , round(rate,2) , " samples/sec    " ) , collapse="" ) )
}

# beta density functions parameterized as prob,theta

dbeta2 <- function( x , prob , theta , log=FALSE ) {
    a <- prob * theta
    b <- (1-prob) * theta
    dbeta( x , shape1=a , shape2=b , log=log )
}

rbeta2 <- function( n , prob , theta ) {
    a <- prob * theta
    b <- (1-prob) * theta
    rbeta( n , shape1=a , shape2=b )
}

# gamma-poisson functions

dgamma2 <- function( x , mu , scale , log=FALSE ) {
    dgamma( x , shape=mu / scale , scale=scale , log=log )
}

rgamma2 <- function( n , mu , scale ) {
    rgamma( n , shape=mu/scale , scale=scale )
}

dgampois <- function( x , mu=NULL , shape=NULL , scale=NULL , rate=NULL , log=FALSE ) {
    if ( !is.null(rate) ) scale <- 1/rate
    if ( !is.null(mu) ) shape <- mu / scale
    prob <- 1 / ( 1 + scale )
    dnbinom( x , size=shape , prob=prob , log=log )
}

rgampois <- function( n , mu=NULL , shape=NULL , scale=NULL , rate=NULL ) {
    if ( !is.null(rate) ) scale <- 1/rate
    if ( !is.null(mu) ) shape <- mu / scale
    prob <- 1 / ( 1 + scale )
    rnbinom( n , size=shape , prob=prob )
}

# mcmc utilities

plotchain <- function(x , mean=TRUE , main=thecall[2][[1]] ) {
    thecall <- match.call()
    plot(x,col="slateblue",type="l",ylab="value",main=main)
    if ( mean==TRUE ) {
        imean <- rep(0,length(x))
        imean[1] <- x[1]
        for ( i in 2:length(x) ) {
            imean[i] <- ( imean[i-1]*(i-1) + x[i] )/i
        }
        lines( 1:length(x) , imean , col="white" , lwd=4 , type="l" )
        lines( 1:length(x) , imean , col="orange" , lwd=2 , type="l" )
    }
}

plotpost <- function(x , burn=0.1 , width=1 , main=thecall[2][[1]] , ci=FALSE , ... ) {
    thecall <- match.call()
    xend <- floor( (length(x)*burn) )
    plot( density( x[-c(1:xend)] , adjust=width ) , col="slateblue" , lwd=2 , main=main , ... )
    if ( ci==TRUE ) {
        x2 <- x[-c(1:xend)] # removed burn in
        mu <- mean( x2 )
        ci95 <- quantile( x2 , probs=c(0.025,0.975) )
        lines( c(mu,mu) , c(0,1000) , col="orange" , lty=2 )
        lines( c(ci95[1],ci95[1]) , c(0,1000) , col="orange" , lty=2 )
        lines( c(ci95[2],ci95[2]) , c(0,1000) , col="orange" , lty=2 )
        themode <- chainmode( x2 )
        lines( c(themode,themode) , c(0,1000) , col="orange" , lty=1 )
    }
}

plotmcmc <- function( mcmc , ci=TRUE , width=1 ) {
    # takes associated chains in the form of a data frame,
    # rows: parameters
    # col 1: chains
    # col 2: posteriors
    npar <- ncol(mcmc)
    par(mfrow=c(npar,2))
    for ( i in 1:npar ) {
        plotchain( mcmc[,i] , main=colnames(mcmc)[i] )
        plotpost( mcmc[,i] , main=colnames(mcmc)[i] , ci=ci , width=width )
    }
}

covmat <- function( m , digits=4 ) {
    # upper diag is covariances
    # lower diag is correlations
    if ( class(m)[1]=="data.frame" ) mcov <- cov( m ) else mcov <- vcov(m)
    mcor <- cov2cor( mcov )
    mcov[ lower.tri(mcov) ] <- NA
    mcor[ lower.tri(mcor) ] <- NA
    result <- list( vcov=round(mcov,digits=digits) , cor=round(mcor,digits=digits) )
    result
}

###
# zero-inflated gamma functions
dzerogamma <- function( x , prob , mu , scale , log=TRUE ) {
    p1 <- dgamma2( x , mu=mu , s=scale , log=FALSE )
    p <- ifelse( x==0 , prob , (1-prob)*p1 )
    if ( log==TRUE ) p <- log(p)
    p
}

####
# drawing functions

arrowspline <- function( from , to , by , shape=c(-1,-1,-1) , arrow=TRUE , arrowlen=0.1 , label=NULL , pos=3 , stem=4 , ... ) {
    if ( class(by)=="matrix" ) {
        # rows for points, col 1 for x, col 2 for y
        xby <- by[,1]
        yby <- by[,2]
    } else {
        xby <- by[1]
        yby <- by[2]
    }
    xs <- xspline( c( from[1] , xby , to[1] ) , c( from[2] , yby , to[2] ) , shape=shape , draw=FALSE )
    lines( xs$x , xs$y , ... )
    n <- length(xs$x)
    arrows( xs$x[n-stem] , xs$y[n-stem] , xs$x[n] , xs$y[n] , length=arrowlen , ... )
    if ( !is.null(label) ) text( xby[1] , yby[1] , label=label , pos=pos , ... )
}

segmentsby <- function( x , y , by , ... ) {
    byid <- unique( by )
    for ( i in byid ) {
        x0 <- x[ by==i ]
        y0 <- y[ by==i ]
        lines( x0 , y0 , ... )
    }
}

shade <- function( object , lim , label=NULL , col="#00000066" , border=NA , ... ) {
    if ( missing(lim) ) stop( "Interval limits missing." )
    if ( missing(object) ) stop( "No density or formula object." )
    from <- lim[1]
    to <- lim[2]
    if ( class(object)=="formula" ) {
        # formula input
        x1 <- eval( object[[3]] )
        y1 <- eval( object[[2]] )
        x <- x1[ x1>=from & x1<=to ]
        y <- y1[ x1>=from & x1<=to ]
    }
    if ( class(object)=="density" ) {
        # density input
        x <- object$x[ object$x>=from & object$x<=to ]
        y <- object$y[ object$x>=from & object$x<=to ]
    }
    polygon( c( x , to , from ) , c( y , 0 , 0 ) , col=col , border=border , ... )
    if ( !is.null(label) ) {
        lx <- mean(x)
        ly <- max(y)/2
        text( lx , ly , label )
    }
}

### softmax rule
multilogistic <- function( x , lambda=1 , diff=TRUE , log=FALSE ) {
    if ( diff==TRUE ) x <- x - min(x)
    f <- exp( lambda * x )
    if ( log==FALSE ) result <- f / sum(f)
    if ( log==TRUE ) result <- log(f) - log(sum(f))
    result
}