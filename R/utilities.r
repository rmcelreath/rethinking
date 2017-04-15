# misc utilities

# various utility functions

# set help to html
#' @export
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


#' Progress display
#'
#' Provides a progress display with estimated time remaining, assuming rate
#' constant process.
#'
#' This function provides useful progress information and estimated time until
#' completion for long looped operations, like sampling from MCMC.
#'
#' @param current Current loop index of process
#' @param min Minimum loop index, usually 0 or 1
#' @param max Maximum loop index
#' @param starttime Time stamp when process began, from \code{Sys.time}
#' @param update.interval How often to display progress
#' @author Richard McElreath
#' @export
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
    mcov <- if ( inherits(m, "data.frame") ) cov(m) else vcov(m)
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


#' Find mode of a continuous density estimate
#'
#' Returns estimated mode of a density computed from samples.
#'
#' This function just finds the x value that maximizes the y density in the
#' density estimate.
#'
#' @param chain Values, e.g. sampled from a posterior via MCMC
#' @param ... Optional arguments passed to density calculation
#' @author Richard McElreath
#' @export
chainmode <- function( chain , ... ) {
    dd <- density(chain , ...)
    dd$x[which.max(dd$y)]
}

# highest posterior density interval, sensu Box and Tiao
# requires coda library
PIprimes <- c(0.67,0.89,0.97) # my snarky prime valued percentiles


#' Confidence/credible intervals from samples
#'
#' These functions compute highest posterior density (HPDI) and percentile (PI)
#' intervals, using samples from a posterior density or simulated outcomes.
#'
#' Highest Posterior Density Intervals (HPDI) are calculated by
#' \code{\link{HPDinterval}} in the \code{coda} package.
#'
#' Percentile intervals (PI) use \code{\link{quantile}} and assign equal mass
#' to each tail.
#'
#' @aliases PI HPDI PCI
#' @param samples Vector of parameter values
#' @param prob interval probability mass
#' @author Richard McElreath
#' @seealso \code{\link{HPDinterval}}
#' @export
HPDI <- function( samples , prob=0.89 ) {
    # require(coda)
    ## class.samples <- class(samples)[1]  # not needed now that we are using inherits()
    coerce.list <- c( "numeric" , "matrix" , "data.frame" , "integer" , "array" )
    if ( inherits(samples, coerce.list) ) {
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
#' @export
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
#' @export
PI <- PCI


#' @export
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
#' @export
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


#' Multi-core version of replicate
#'
#' Uses the \code{parallel} library to distribute \code{\link{replicate}}
#' processing across cores.
#'
#' This function uses \code{\link{mclapply}} to distribute replications across
#' cores. It then simplifies the result to an array.
#'
#' @param n Number of replications
#' @param expr Code to replicate
#' @param refresh Status update refresh interval
#' @param mc.cores Number of cores to use
#' @author Richard McElreath
#' @seealso \code{\link{mclapply}}, \code{\link{replicate}}
#' @export
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
#' @export
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


#' Build and check integer index variables
#'
#' These functions assist with building (\code{coerce_index}) and checking
#' (\code{check_index}) integer index variables of the kind needed to define
#' varying effect models.
#'
#' Varying effect models often require index variables that begin at 1 and
#' comprise only integers. These variables are used to lookup specific
#' parameter values inside of the model. For example, it is common to define
#' varying intercepts with a linear model like \code{a0 + a_id[id[i]]}. Here
#' the variable \code{id} is an index variable. It has one value for each case,
#' defining which individual applies to that case.
#'
#' When raw data exist as factors, these index variables much be converted to
#' integers. This is trickier than it sounds, because R uses an internal
#' integer represntation for factors, \code{levels}, that can conflict with
#' ordinary integer representations.
#'
#' The function \code{coerce_index} deals with that complication. When the
#' input is a single vector of factors, it returns an integer vector beginning
#' at 1 and with contiguous values.
#'
#' When the input is instead a comma-separated list of factors, it returns a
#' list in which each factor has been converted to integers, but all levels in
#' all factors were merged so that the same labels across factors always have
#' the same integer values in the result. For example, suppose cases refer to
#' dyads and there are two factors, \code{id1} and \code{id2}, that indicate
#' which pair of individuals are present in each dyad. The labels in these
#' variables should refer to the same individuals. Passing both simultaneously
#' to \code{coerce_index} ensures that the results respect that fact.
#'
#' The function \code{check_index} merely checks an integer vector to see if it
#' is contiguous.
#'
#' @aliases coerce_index check_index
#' @param ... A comma-separted list of variables. See details.
#' @param x A vector of integers to check for contiguity
#' @return For \code{coerce_index}, the result is either a single vector of
#' integers (if the input was a single vector) or rather a list of vectors of
#' integers (if the input was a list of vectors).
#' @author Richard McElreath
#' @export
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
