\name{cherry_blossoms}
\alias{cherry_blossoms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Japan Cherry Blossom Historical Data}
\description{
  Historical Series of Phenological data for Cherry Tree Flowering at Kyoto City.
}
\usage{
data(cherry_blossoms)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
}
\format{
    \enumerate{
        \item year: Year CE
        \item doy: Day of year of first bloom. Day 89 is April 1. Day 119 is May 1.
        \item temp: March temperature estimate
        \item temp_upper: Upper 95\% bound for estimate
        \item temp_lower: Lower 95\% bound for estimate 
    }
}
\value{
}
\references{
Aono and Saito 2010. International Journal of Biometeorology, 54, 211-219.
Aono and Kazui 2008. International Journal of Climatology, 28, 905-914.
Aono 2012. Chikyu Kankyo (Global Environment), 17, 21-29.
http://atmenv.envi.osakafu-u.ac.jp/aono/kyophenotemp4/}
\seealso{}
\examples{

# This code reproduces the plot on the 2nd edition cover

library(rethinking)
data(cherry_blossoms)
d <- cherry_blossoms

# spline on temp
d2 <- d[ complete.cases(d$temp) , ] # complete cases on temp

num_knots <- 30
( knot_list <- quantile( d2$year , probs=seq(0,1,length.out=num_knots) ) )

library(splines)
B <- bs(d2$year, 
    knots=knot_list[-c(1,num_knots)] , 
    degree=3 , intercept=TRUE )

m1 <- quap(
    alist(
        T ~ dnorm( mu , sigma ) ,
        mu <- a + B %*% w ,
        a ~ dnorm(6,1),
        w ~ dnorm(0,10),
        sigma ~ dexp(1)
    ),
    data=list( T=d2$temp , B=B ) , 
    start=list( w=rep( 0 , ncol(B) ) ) )

# now spline on blossom doy
d3 <- d[ complete.cases(d$doy) , ] # complete cases on doy

knot_list <- seq( from=min(d3$year) , to=max(d3$year) , length.out=num_knots )
B3 <- t(bs(d3$year, knots=knot_list , degree=3, intercept = FALSE))

m2 <- quap(
    alist(
        Y ~ dnorm( mu , sigma ) ,
        mu <- a0 + as.vector( a %*% B ),
        a0 ~ dnorm(100,10),
        a ~ dnorm(0,10),
        sigma ~ dexp(1)
    ),
    data=list( Y=d3$doy , B=B3 ) , start=list(a=rep(0,nrow(B3))) )

# PLOT

blank2(w=2,h=2)

par( mfrow=c(2,1) , mgp = c(1.25, 0.25, 0), mar = c(0.75, 2.5, 0.75, 0.75) + 0.1, 
        tck = -0.02, cex.axis = 0.8 )

xcex <- 1.2
xpch <- 16
xcol1 <- col.alpha(rangi2,0.3)
col_spline <- col.alpha("black",0.4)
xlims <- c(850,2000)

plot( d2$year , d2$temp , ylab="March temperature" , col=xcol1 , pch=xpch , cex=xcex , xlab="" , xlim=xlims , bty="n" , axes=FALSE , ylim=c( 4.5 , 8.3 ) )
l <- link( m1 )
li <- apply(l,2,PI,0.97)

atx <- c(900,1400,2000)
axis( 1 , at=atx , labels=paste(atx,"CE") )
axis( 2 , at=c(5,8) , labels=c("5°C","8°C") )

y <- d3$doy
y <- y - min(y)
y <- y/max(y)
blossom_col <- sapply( d3$doy , function(y) hsv(1, rbeta2(1, inv_logit(logit(0.1)+0.02*y) ,10) ,1,0.8) )
plot( NULL , cex=xcex , ylab="Day of first blossom" , xlim=xlims , bty="n" , axes=FALSE , xlab="" , ylim=range(d3$doy) )
l <- link( m2 )
li <- apply(l,2,PI,0.9)
points( d3$year , d3$doy , col=blossom_col ,  pch=8  , cex=xcex , lwd=2 )
shade( li , d3$year , col=grau(0.3) )

axis( 2 , at=c(90,120) , labels=c("April 1","May 1") )

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

