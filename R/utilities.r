# misc utilities

# various utility functions

# default quartz plot size for book: 3.5in by 4in, giving square plot for default margins
blank <- function(ex=1) quartz("myquartz",width=3.5*ex,height=4*ex)
blank2 <- function() {
    blank(ex=0.9)
    quartzFonts(serif = quartzFont(rep("MinionPro-Regular", 4)))
    par(family="serif")
}

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

Rho <- function( model , digits=2 ) {
    round( cov2cor(vcov(model)) , digits )
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

# finds mode of a continuous density
chainmode <- function( chain , ... ) {
    dd <- density(chain , ...)
    dd$x[which.max(dd$y)]
}

# highest posterior density interval, sensu Box and Tiao
# requires coda library
HPDI <- function( samples , prob=0.95 ) {
    # require(coda)
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

se <- function( model ) {
    sqrt( diag( vcov(model) ) )
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

# replicate with progress display
replicate2 <- function (n, expr, interval=0.1, simplify = "array") {
    show_progress <- function(i) {
        intervaln <- floor( n * interval )
        if ( floor(i/intervaln) == i/intervaln ) {
            cat( paste( "[" , i , "/" , n , "]\r" ) )
        }
    }
    result <- sapply(1:n, 
        eval.parent(substitute(function(i,...) { show_progress(i); expr })), 
        simplify = simplify)
    cat("\n")
    result
}
