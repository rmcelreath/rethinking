# plotting

#' @export
set_nice_margins <- function() {
    par_mf <- par("mfrow","mfcol")
    if ( all(unlist(par_mf)==1) ) {
        #par_old <- par(no.readonly = TRUE)
        #on.exit(par(par_old))
        par(mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1, tck = -0.02)
    }
}



#' Density plots
#'
#' Convenient interface for plotting density estimates.
#'
#' This function merely provides a convenient interface for plotting density
#' estimates produced by \code{density}. It handles both single vectors and
#' multiple vectors, contained within a data frame.
#'
#' Highest Posterior Density Intervals (HPDI) are calculated by
#' \code{HPDinterval} in the \code{coda} package.
#'
#' @param x Vector of values to construct density from, or data frame. If
#' \code{x} is a data frame, then \code{dens} plots a grid of densities, one
#' for each column in \code{x}.
#' @param adj width of density kernal.
#' @param norm.comp If \code{TRUE}, overlays normal density comparison.
#' @param show.HPDI If a numeric value, displays HPDI of same width. For
#' example, \code{show.HPDI=0.95} shows a 95 percent HPDI inside the density.
#' @param show.zero If \code{TRUE}, draws a vertical line at location of zero
#' on horizonal axis.
#' @param rm.na If \code{TRUE}, removes \code{NA}s before computing density
#' @param add If \code{TRUE}, overlays density on an existing plot
#' @param ... Other parameters to pass to \code{density}, which constructs the
#' density estimates.
#' @author Richard McElreath
#' @seealso \code{\link{density}}, \code{\link{HPDinterval}}
#' @export
#'
dens <- function( x , adj=0.5 , norm.comp=FALSE , main="" , show.HPDI=FALSE , show.zero=FALSE , rm.na=TRUE , add=FALSE , ...) {
    if ( inherits(x, "data.frame") ) {
        # full posterior
        n <- ncol(x)
        cnames <- colnames(x)
        set_nice_margins()
        par( mfrow=make.grid(n) )
        for ( i in 1:n ) {
            dens( x[,i] , adj=adj , norm.comp=norm.comp , show.HPDI=show.HPDI , show.zero=TRUE , xlab=cnames[i] , ... )
        }
    } else {
        # vector
        if ( rm.na==TRUE ) x <- x[ !is.na(x) ]
        thed <- density(x,adjust=adj)
        if ( add==FALSE ) {
            set_nice_margins()
            plot( thed , main=main , ... )
        } else
            lines( thed$x , thed$y , ... )
        if ( show.HPDI != FALSE ) {
            hpd <- HPDI( x , prob=show.HPDI )
            shade( thed , hpd )
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


#' Contour plot from equal length x,y,z vectors
#'
#' Provides an interface to use \code{contour} by providing three equal length
#' vectors for x, y and z coordinates.
#'
#' This function merely constructs a matrix suitable for \code{contour}, using
#' x, y and z coordinates.
#'
#' @param x vector of x values
#' @param y vector of y values
#' @param z vector of z values
#' @param ... other parameters to pass to \code{contour}
#' @author Richard McElreath
#' @seealso \code{\link{contour}}
#' @export
#'
contour_xyz <- function( x , y , z , ... ) {
    ux <- unique(x)
    uy <- unique(y)
    n <- length(ux)
    m <- matrix( z , nrow=n , ncol=n )
    contour( ux , uy , m , ... )
}

# just converts inputs to form expected by image()


#' Heat map from equal length x,y,z vectors
#'
#' Provides an interface to use \code{image} by providing three equal length
#' vectors for x, y and z coordinates.
#'
#' This function merely constructs a matrix suitable for \code{image}, using x,
#' y and z coordinates.
#'
#' @param x vector of x values
#' @param y vector of y values
#' @param z vector of z values
#' @param ... other parameters to pass to \code{image}
#' @author Richard McElreath
#' @seealso \code{\link{image}}
#' @export
#'
image_xyz <- function( x , y , z , ... ) {
    image( unique(x) , unique(y) , matrix(z, length(unique(x)), length(unique(y)) ) , ... )
}

# plot new y values on secondary y-axis on right side
#' @export
plot2y <- function( ... , y2lab="2nd axis" , y2col=NULL ) {
    if ( is.null(y2col) ) y2col <- "black"
    par(new=TRUE)
    par(mar=c(5,4,4,5)+.1)
    plot( xaxt="n" , yaxt="n" , xlab="" , ylab="" , ... )
    axis( 4 , col=y2col , col.axis=y2col )
    mtext( y2lab , side=4 , line=3 , col=y2col )
}


# function for plotting "naive" posteriors and 95% confints
#' @export
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


#' Simple histograms
#'
#' Simple integer-valued histograms, for displaying count distributions.
#'
#' This function constructs clean histograms for count data. Non-integer data
#' can be rounded to nearest integer before plotting. Internally, this function
#' is little more than \code{plot(table(x))}.
#'
#' @param x Vector of values to construct histogram from
#' @param round When \code{TRUE}, rounds values in \code{x} before plotting
#' @param ylab Label on vertical axis
#' @param ... Other parameters to pass to plot
#' @author Richard McElreath
#' @seealso \code{\link{hist}}
#' @export
#'
simplehist <- function( x , round=TRUE , ylab="Frequency" , ... ) {
    if ( round==TRUE ) x <- round(x)
    plot(table(x),ylab=ylab,...)
}

simplehist_old <- function( x , ylab="Frequency" , xlab="Count" , ycounts=TRUE , adjust=1 , lcol="black" , bins=NULL , show.counts=0 , xlim=NULL , ylim=NULL , ... ) {
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
        set_nice_margins()
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
       set_nice_margins()
       plot( density( x , adjust=adjust ) , ylab=ylab , xlab=xlab , ... )
    }

}

####
# drawing functions

#' @export
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

#' @export
segmentsby <- function( x , y , by , ... ) {
    byid <- unique( by )
    for ( i in byid ) {
        x0 <- x[ by==i ]
        y0 <- y[ by==i ]
        lines( x0 , y0 , ... )
    }
}



#' Shade density intervals
#'
#' Adds a shaded interval to an existing density plot or regression plot.
#'
#' This function uses \code{\link{polygon}} to draw a shaded region under a
#' density curve or on a regression plot. The function assumes the plot already
#' exists. See the examples below.
#'
#' There are two plotting interfaces, for densities. First, if the density is
#' plotted from kernal estimation, using perhaps \code{\link{density}}, then
#' the same density estimate should be passed as the first parameter. See
#' example. Second, if the density is plotted from a series of x and y values,
#' from perhaps a grid approximation, then a formula can be passed to define
#' the curve. See example.
#'
#' For plotting confidence regions on a regression plot, the matrix object
#' should contain y-axis values defining the border of the region. The first
#' row of the matrix should be the bottom of the region, and the second row
#' should be the top. Each column should correspond to the x-axis values in
#' \code{lim}. See example.
#'
#' @param object A \code{density} or \code{formula} object that defines the
#' density OR a \code{matrix} object that defines the plot region. See details.
#' @param lim For a density, a vector of two values, indicating upper and lower
#' bounds of interval, on the horizontal axis. For a plot region, a list of
#' x-axis values corresponding to y-axis points in the matrix object.
#' @param label Optional label to center in interval.
#' @param col Color to shade the interval. Default is transparent gray.
#' @param border Border of shaded region. Default is no border, \code{NA}.
#' @param ... Other parameters to pass to \code{polygon}, which actually draws
#' the shaded region.
#' @author Richard McElreath
#' @seealso \code{\link{density}}, \code{\link{dens}}
#' @export
#' @examples
#'
#' models <- seq( from=0 , to=1 , length.out=100 )
#' prior <- rep( 1 , 100 )
#' likelihood <- dbinom( 6 , size=9 , prob=models )
#' posterior <- likelihood * prior
#' posterior <- posterior / sum(posterior)
#'
#' # using formula interface
#' plot( posterior ~ models , type="l" )
#' shade( posterior ~ models , c(0,0.5) )
#'
#' # using density interface
#' samples <- sample( models , size=1e4 , replace=TRUE , prob=posterior )
#' plot( density(samples) )
#' shade( density(samples) , PCI(samples) )
#'
#' # plotting a shaded confidence interval on a regression plot
#' data(cars)
#' m <- lm( dist ~ speed , cars )
#' p <- extract.samples( m )
#' x.seq <- seq( from=min(cars$speed)-1 , to=max(cars$speed)+1 , length.out=30 )
#' mu.ci <- sapply( x.seq , function(x) PI( p[,1] + p[,2]*x ) )
#' plot( dist ~ speed , cars )
#' abline( m )
#' shade( mu.ci , x.seq )
#'
shade <- function( object , lim , label=NULL , col=col.alpha("black",0.15) , border=NA , ... ) {
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
    if ( class(object)=="matrix" & length(dim(object))==2 ) {
        # matrix defining confidence region around a curve
        y <- c( object[1,] , object[2,][ncol(object):1] ) # reverse second row
        x <- c( lim , lim[length(lim):1] ) # lim needs to be x-axis values
    }
    # draw
    if ( class(object)=="matrix" ) {
        polygon( x , y , col=col , border=border , ... )
    } else {
        polygon( c( x , to , from ) , c( y , 0 , 0 ) , col=col , border=border , ... )
    }
    # label?
    if ( !is.null(label) ) {
        lx <- mean(x)
        ly <- max(y)/2
        text( lx , ly , label )
    }
}

#' @export
mcmcpairs <- function( posterior , cex=0.3 , pch=16 , col=col.alpha("slateblue",0.2) , n=1000 , adj=1 , ... ) {
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
        points( x[i] , y[i] , ... )
    }
    panel.cor <- function( x , y , ... ) {
        k <- cor( x , y )
        cx <- sum(range(x))/2
        cy <- sum(range(y))/2
        text( cx , cy , round(k,2) , cex=2*exp(abs(k))/exp(1) )
    }
    set_nice_margins()
    pairs( posterior , cex=cex , pch=pch , col=col , upper.panel=panel.2d , lower.panel=panel.cor , diag.panel=panel.dens , ... )
}

