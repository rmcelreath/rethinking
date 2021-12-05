# my DAG drawing function to extend dagitty plot method

drawdag <- function( x , col_arrow="black" , col_segment="black" , col_labels="black" , cex=1 , lwd=1.5 , goodarrow=TRUE , xlim , ylim , shapes , col_shapes , radius=3.5 , add=FALSE , xkcd=FALSE , latent_mark="c" , ... ){ 
    require(dagitty)

    # check for list of DAGs
    if ( class(x)=="list" ) {
        n <- length(x)
        y <- make.grid(n)
        par(mfrow=y)
        for ( i in 1:n ) drawdag( x[[i]] , ... )
        return(invisible(NULL))
    }

    x <- as.dagitty( x )
    dagitty:::.supportsTypes(x,c("dag","mag","pdag"))
    coords <- coordinates( x )
    if( any( !is.finite( coords$x ) | !is.finite( coords$y ) ) ){
        coords <- coordinates( graphLayout( x ) )
        #stop("Please supply plot coordinates for graph! See ?coordinates and ?graphLayout.")
    }
    labels <- names(coords$x)
    if ( add==FALSE ) {
        par(mar=rep(0,4))
        plot.new()
        par(new=TRUE)
    }
    wx <- sapply( paste0("mm",labels), 
        function(s) strwidth(s,units="inches") )
    wy <- sapply( paste0("\n",labels), 
        function(s) strheight(s,units="inches") )
    ppi.x <- dev.size("in")[1] / (max(coords$x)-min(coords$x))
    ppi.y <- dev.size("in")[2] / (max(coords$y)-min(coords$y))
    wx <- wx/ppi.x
    wy <- wy/ppi.y
    if ( missing(xlim) )
        xlim <- c(min(coords$x-wx/2),max(coords$x+wx/2))
    if ( missing(ylim) )
        ylim <- c(-max(coords$y+wy/2),-min(coords$y-wy/2))
    if (add==FALSE ) plot( NA, xlim=xlim, ylim=ylim, xlab="", ylab="", bty="n",
        xaxt="n", yaxt="n" )
    # for each variable marked latent, set to draw with circle
    unobs_vars <- latents(x)
    if ( length(unobs_vars)>0 & !is.null(latent_mark) ) {
        if ( missing(shapes) ) shapes <- list()
        for ( uv in unobs_vars ) {
            if ( is.null(shapes[[uv]]) )
                shapes[[ uv ]] <- latent_mark
        }
    }
    # review node labels for subscript _
    xlabels <- lapply( labels , function(l) {
            lr <- l
            xx <- strsplit( l , "_" )[[1]]
            if ( length(xx) > 1 ) {
                lr <- bquote( .(xx[1])[.(xx[2])] )
            }
            return(lr)
        } )
    xxbuffer <- strwidth("xx")
    yybuffer <- strheight("\n")
    wx <- sapply( xlabels, 
        function(s) strwidth(s) + xxbuffer )
    wy <- sapply( xlabels,
        function(s) max( strheight(s) , yybuffer ) )
    names(wx) <- labels
    names(wy) <- labels
    if ( length(unobs_vars)>0 & !is.null(latent_mark) ) {
        # check for latent vars and increase space buffer
        wx[ unobs_vars ] <- xxbuffer*2.1
        wy[ unobs_vars ] <- yybuffer*1.1
    }
    # edges
    asp <- par("pin")[1]/diff(par("usr")[1:2]) /
        (par("pin")[2]/diff(par("usr")[3:4]))
    ex <- dagitty::edges(x)
    ax1 <- rep(0,nrow(ex))
    ax2 <- rep(0,nrow(ex))
    ay1 <- rep(0,nrow(ex))
    ay2 <- rep(0,nrow(ex))
    axc <- rep(0,nrow(ex))
    ayc <- rep(0,nrow(ex))
    acode <- rep(2,nrow(ex))
    has.control.point <- rep(FALSE,nrow(ex))
    for( i in seq_len(nrow(ex)) ){
        if( ex[i,3] == "<->" ){
            acode[i] <- 3
            has.control.point[i] <- TRUE
        }
        if( ex[i,3] == "--" ){
            acode[i] <- 0
        }
        l1 <- as.character(ex[i,1]); l2 <- as.character(ex[i,2])
        x1 <- coords$x[l1]; y1 <- coords$y[l1]
        x2 <- coords$x[l2]; y2 <- coords$y[l2]
        if( is.na( ex[i,4] ) || is.na( ex[i,5] ) ){
            cp <- dagitty:::.autoControlPoint( x1, y1, x2, y2, asp,
                .2*as.integer( acode[i]==3 ) )
        } else {
            cp <- list(x=ex[i,4],y=ex[i,5])
            has.control.point[i] <- TRUE
        }
        bi1 <- dagitty:::.lineSegBoxIntersect( x1-wx[l1]/2,y1-wy[l1]/2,
            x1+wx[l1]/2,y1+wy[l1]/2, x1, y1, cp$x, cp$y )
        bi2 <- dagitty:::.lineSegBoxIntersect( x2-wx[l2]/2,y2-wy[l2]/2,
            x2+wx[l2]/2,y2+wy[l2]/2, cp$x, cp$y, x2, y2 )
        if( length(bi1) == 2 ){
            x1 <- bi1$x; y1 <- bi1$y
        }
        if( length(bi2) == 2 ){
            x2 <- bi2$x; y2 <- bi2$y
        }
        ax1[i] <- x1; ax2[i] <- x2
        ay1[i] <- y1; ay2[i] <- y2
        axc[i] <- cp$x; ayc[i] <- cp$y
    }
    directed <- acode==2 & !has.control.point
    undirected <- acode==0 & !has.control.point

    arr.width <- 0.15
    arr.type <- "curved"
    arr.adj <- 1
    if ( xkcd==TRUE ) {
        for ( ii in 1:length(ax1[directed]) ) {
            lines_xkcd( c(ax1[directed][ii],ax2[directed][ii]) , c(-ay1[directed][ii],-ay2[directed][ii]) , col=col_segment , lwd=lwd*2 , lwdbg=lwd*4 , seg=10 )
        }#ii
        arr.width <- arr.width * lwd
        arr.type <- "triangle"
        arr.adj <- 0.5
        goodarrow <- TRUE
    }

    if ( goodarrow==TRUE ) {
        #require(shape)
        shape::Arrows( ax1[directed], -ay1[directed], 
            ax2[directed], -ay2[directed], arr.length=0.2 , arr.width=arr.width, col=col_arrow , lwd=lwd , arr.adj=arr.adj , arr.type=arr.type )
    } else
        arrows( ax1[directed], -ay1[directed], 
            ax2[directed], -ay2[directed], length=0.1, col=col_arrow , lwd=lwd )

    segments( ax1[undirected], -ay1[undirected], ax2[undirected], -ay2[undirected], col=col_segment , lwd=lwd )
    
    for( i in which( has.control.point ) ){
        dag_arc( ax1[i], -ay1[i], 
            ax2[i], -ay2[i], axc[i], -ayc[i], 
            col=c( col_arrow , col_segment )[1+(acode[i]==0)], 
            code=acode[i], length=0.1, lwd=lwd+(acode[i]==0) , goodarrow=goodarrow )
    }
    
    # node shapes?
    # should be named list with "c" for circle or "b" for box
    if ( !missing(shapes) ) {
        if ( length(shapes)>0 ) {
        for ( i in 1:length(shapes) ) {
            the_label <- names(shapes)[i]
            j <- which( labels==the_label )
            if ( missing(col_shapes) ) col_shapes <- col_labels[ min(j,length(col_labels)) ]
            if ( length(j)>0 ) {
                cpch <- 1
                if ( shapes[[i]]=="fc" ) cpch <- 16
                if ( shapes[[i]] %in% c("c","fc") ) 
                    #circle( coords$x[the_label] , -coords$y[the_label] , r=radius , lwd=lwd , col=col_shapes )
                    # draw circle with point so it has correct aspect ratio
                    points( coords$x[the_label] , -coords$y[the_label] , cex=radius , lwd=lwd , col=col_shapes , pch=cpch )
            }
        }#i
        }#>0
    }
    # node labels
    for ( i in 1:length(xlabels) ) {
        k <- col_labels
        #if ( length(col_labels)>1 ) k <- col_labels[i]
        text( coords$x[i] , -coords$y[labels][i] , xlabels[[i]] , cex=cex , col=k )
    }

    # prep arrow coordinate table
    arrtab <- data.frame( ax1 , ay1 , ax2 , ay2 , directed )
    return( invisible(arrtab) )
}

circle <- function( x , y , r=1 , npts=100 , ... ) {
    theta <- seq( 0, 2*pi , length = npts )
    lines( x = x + r * cos(theta) , y = y + r * sin(theta) , ... )
}

dag_arc <- function (x1, y1, x2, y2, xm, ym, col = "gray", length = 0.1, 
    code = 3, lwd = 1.5 , goodarrow=TRUE ) {
    x <- c(x1, xm, x2)
    y <- c(y1, ym, y2)
    res <- xspline(x, y, 1, draw = FALSE)
    lines( res , col = col , lwd=lwd )
    nr <- length(res$x)
    if (code >= 3) {
        if ( goodarrow==TRUE ) {
            #require(shape)
            shape::Arrows( res$x[4], res$y[4], res$x[1], res$y[1], arr.length=0.2 , arr.width=0.15, col=col , lwd=lwd , arr.adj=1 , arr.type="curved" )
        } else
            arrows(res$x[1], res$y[1], res$x[4], res$y[4], col = col, 
            code = 1, length = length, lwd = lwd)
    }
    if (code >= 2) {
        if ( goodarrow==TRUE ) {
            #require(shape)
            shape::Arrows( res$x[nr-3], res$y[nr-3], res$x[nr], res$y[nr], arr.length=0.2 , arr.width=0.15, col=col , lwd=lwd , arr.adj=1 , arr.type="curved" )
        } else
            arrows(res$x[nr - 3], res$y[nr - 3], res$x[nr], res$y[nr], 
            col = col, code = 2, length = length, lwd = lwd)
    }
}

# function to map coordinates from dag x to reduced dag y
# returns new dagitty
# used by drawopenpaths
dag_copy_coords <- function( x , y ) {
    orig <- coordinates(x)
    in_dag <- names( coordinates(y)$x )
    new_coords <- coordinates(y)
    new_coords$x <- orig$x[ in_dag ]
    new_coords$y <- orig$y[ in_dag ]
    coordinates(y) <- new_coords
    return(y)
}

# function to overlay open paths, conditioning on list Z
drawopenpaths <- function( x , Z=list() , col_arrow="red" , ... ) {
    x <- as.dagitty( x )
    # get all paths
    path_list <- paths( x , Z=Z )
    shapes <- list()
    if ( length(Z)>0 ) {
        for ( i in 1:length(Z) ) shapes[[ Z[[i]] ]] <- "c"
    }
    # draw open paths in highlight color
    for ( i in 1:length(path_list) ) {
        if ( path_list$open[i]==TRUE ) {
            path_dag <- dagitty( concat( "dag { " , path_list$paths[i] , " }" ) )
            path_dag <- dag_copy_coords( x , path_dag )
            drawdag( path_dag , col_arrow=col_arrow , add=TRUE , ... )
        }
    }#i
    return(invisible(path_list))
}

# test code
if (FALSE) {

plot( NULL , xlim=c(-1,1) , ylim=c(-1,1) )
circle( 0 , 0 , 1 )

library(rethinking)
library(dagitty)
plant_dag <- dagitty( "dag {
    h0 -> h1
    f -> h1
    t -> f
}")
coordinates( plant_dag ) <- list( x=c(h0=0,t=2,f=1,h1=1) , y=c(h0=0,t=0,f=1,h1=2) )

drawdag( plant_dag , cex=1.2 , col_labels=c("red","black","red","red") , col_arrow=c("red","black","red") , goodarrow=TRUE )

exdag <- dagitty( "dag {
    U [unobserved]
    Z -> X -> Y
    X <- U -> Y
}")
coordinates( exdag ) <- list( x=c(Z=0,X=1,Y=2,U=1.5) , y=c(Z=0,X=0,Y=0,U=-1) )
drawdag( exdag , radius=3.8 )

# drawing paths

g <- dagitty( "dag { 
    x -> y
    a -> x
    b -> y
    a -> z <- b
    x [exposure]
    y [outcome] 
}" , layout=TRUE )

drawdag( g , col_arrow="gray" )

drawopenpaths( g , col_arrow="black" )
drawopenpaths( g , Z=list("z") , col_arrow="black" )

# xkcd example
exdag <- dagitty( "dag {
    z -> x -> y
    x <- U -> y
}")
coordinates( exdag ) <- list( x=c(z=0,x=1,y=2,U=1.5) , y=c(z=0,x=0,y=0,U=-1) )
drawdag( exdag , xkcd=TRUE , lwd=1.5 )



}

