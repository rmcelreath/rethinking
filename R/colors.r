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

rangi2 <- col.desat( "blue" , 0.5 )
grau <- function( alpha=0.5 ) col.alpha( "black" , alpha )

#rethink_palette <- c("#5BBCD6","#F98400","#F2AD00","#00A08A","#FF0000")
rethink_palette <- c("#8080FF","#F98400","#00A08A","#E2AD00","#FF0000")
rethink_cmyk <- c(col.alpha("black",0.25),"cyan")
