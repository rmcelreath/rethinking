# misc utilities

# various utility functions

# set help to html
htmlhelp <- function() options(help_type="html")

# set CRAN mirror
setcran <- function(themirror="http://cran.stat.ucla.edu/") options(repos=structure(c(CRAN=themirror)))

# default quartz plot size for book: 3.5in by 4in, giving square plot for default margins
blank <- function(ex=1,w=1,h=1) {
    quartz("myquartz",width=3.5*ex*w,height=3.5*ex*h)
    par(mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1, tck = -0.02)
}

# default pdf plot size, for making cmyk figures
# close file with dev.off() as usual
pdfblank <- function (ex = 1, w = 1, h = 1, colormodel="cmyk" , ... ) 
{
    pdf("mypdf.pdf", width = 3.5 * ex * w, height = 3.5 * ex * 
        h , colormodel=colormodel , ...)
    par(mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1, 
        tck = -0.02)
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
#glmm <- function( ... , family , REML=FALSE ) {
#    require(lme4)
#    if ( missing(family) ) {
#        result <- lmer( ... , REML=REML )
#    } else {
#        result <- glmer( ... , family=family )
#    }
#    result
#}

#bglmm <- function( ... , REML=FALSE ) {
#    require(lme4)
#    require(blme)
#    blmer( ... , REML=REML )
#}

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

# finds mode of a continuous density
chainmode <- function( chain , ... ) {
    dd <- density(chain , ...)
    dd$x[which.max(dd$y)]
}

# highest posterior density interval, sensu Box and Tiao
# requires coda library
PIprimes <- c(0.67,0.89,0.97) # my snarky prime valued percentiles
HPDI <- function( samples , prob=0.89 ) {
    # require(coda)
    class.samples <- class(samples)[1]
    coerce.list <- c( "numeric" , "matrix" , "data.frame" , "integer" , "array" )
    if ( class.samples %in% coerce.list ) {
        # single chain for single variable
        samples <- coda::as.mcmc( samples )
    }
    x <- sapply( prob , function(p) coda::HPDinterval( samples , prob=p ) )
    # now order inside-out in pairs
    n <- length(prob)
    result <- rep(0,n*2)
    for ( i in 1:n ) {
        low_idx <- n+1-i
        up_idx <- n+i
        # lower
        result[low_idx] <- x[1,i]
        # upper
        result[up_idx] <- x[2,i]
        # add names
        names(result)[low_idx] <- concat("|",prob[i])
        names(result)[up_idx] <- concat(prob[i],"|")
    }
    return(result)
}

# percentile confidence/credible interval
PCI <- function( samples , prob=0.89 ) {
    x <- sapply( prob , function(p) {
        a <- (1-p)/2
        quantile( samples , probs=c(a,1-a) )
    } )
    # now order inside-out in pairs
    n <- length(prob)
    result <- rep(0,n*2)
    for ( i in 1:n ) {
        low_idx <- n+1-i
        up_idx <- n+i
        # lower
        result[low_idx] <- x[1,i]
        # upper
        result[up_idx] <- x[2,i]
        # add names
        a <- (1-prob[i])/2
        names(result)[low_idx] <- concat(round(a*100,0),"%")
        names(result)[up_idx] <- concat(round((1-a)*100,0),"%")
    }
    return(result)
}
PI <- PCI


se <- function( model ) {
    sqrt( diag( vcov(model) ) )
}

# quadratic estimate confidence intervals from means and standard errors
confint_quad <- function( model=NULL , est , se , prob=0.89 ) {
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
    p <- (1-prob)/2
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

# multi-core replicate
mcreplicate <- function (n, expr, refresh = 0.1, mc.cores=2 ) {
    #require(parallel)
    show_progress <- function(i) {
        intervaln <- floor(n * refresh)
        if (floor(i/intervaln) == i/intervaln) {
            cat(paste("[", i, "/", n, "]\r"))
        }
    }
    result <- simplify2array(mclapply(1:n, eval.parent(substitute(function(i, 
        ...) {
        if (refresh>0) show_progress(i)
        expr
    })),mc.cores=mc.cores))
    if (refresh>0) cat("\n")
    result
}

# check index vector
check_index <- function( x ) {
    y <- sort(unique(x))
    n <- length(y)
    message( concat( "Length: ",n ) )
    message( concat( "Range: ",min(y)," / ",max(y) ) )
    if ( max(y) != n ) message( "Maximum index different than number of unique values" )
    diffs <- sapply( 2:n , function(i) y[i] - y[i-1] )
    if ( any(diffs)!=1 ) message( "At least one gap in consecutive values" )
}

# old vector-only coerce_index
# coerce_index <- function( x ) as.integer(as.factor(as.character(x)))

# new coerce_index that can take multiple vectors
# ensures labels in all are part of same index set
coerce_index <- function( ... ) {
    L <- list(...)
    if ( is.list(L[[1]]) && length(L)==1 ) L <- L[[1]]
    if ( length(L)==1 ) {
        return( as.integer(as.factor(as.character(L[[1]]))) )
    } else {
        # multiple inputs
        vnames <- match.call()
        vnames <- as.character(vnames)[2:(length(L)+1)]
        # generate levels that include all labels from all inputs
        M <- L
        for ( i in 1:length(L) ) M[[i]] <- as.character(L[[i]])
        Mall <- M[[1]]
        for ( i in 2:length(L) ) Mall <- c( Mall , M[[i]] )
        Mall <- unique(Mall)
        new_levels <- levels(as.factor(Mall))
        for ( i in 1:length(L) ) {
            M[[i]] <- factor(M[[i]],levels=new_levels)
            M[[i]] <- as.integer(M[[i]])
        }
        names(M) <- paste( vnames , "_idx" , sep="" )
        return(M)
    } 
}

# sd and var functions that don't use n-1 denominator
sd2 <- function( x , na.rm=TRUE ) {
    sqrt( var2(x,na.rm=na.rm) )
}

var2 <- function( x , na.rm=TRUE ) {
    # use E(x^2) - E(x)^2 form
    mean(x^2) - mean(x)^2
}

# function for formatting summary table output in various show methods
# x should be a list or data frame
# digits should be a named vector with digits to display for each column
# name 'default__' is default digits
# so can pass e.g. digits=c( 'default__'=2 , n_eff=0 )
format_show <- function( x , digits ) {
    r <- as.data.frame(lapply( 1:length(x) , 
        function(i) { 
            if ( names(x)[i] %in% names(digits) ) 
                round( x[[i]] , digits[names(x)[i]] ) 
            else 
                round(x[[i]], digits['default__'] );
        } ) )
    names(r) <- names(x)
    rownames(r) <- rownames(x)
    return(r)
}

# just for recoding "Yes"/"No" to 1/0
yn2bin <- function(x) ifelse(is.na(x),NA,ifelse(x %in% c("Yes","yes","Y","y"),1,0))
