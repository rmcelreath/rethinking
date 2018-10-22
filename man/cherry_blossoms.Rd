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

library(rethinking)
data(cherry_blossoms)
d <- cherry_blossoms

par(mfrow=c(2,1))

plot( d$year , d$blossoms_doy , xlab="Year (CE)" , ylab="Date (day-in-year)" , pch=16 , cex=1.5 , col=col.alpha("pink",0.8) , bty="l" , xlim=c(800,2015) , cex.axis=0.8 )
abline( h=91 , lty=2 , lwd=0.5 )
text( 810 , 89 , "1 April" , cex=0.8 )
abline( h=121 , lty=2 , lwd=0.5 )
text( 810 , 119 , "1 May" , cex=0.8 )
mtext( "Date of first bloom" , adj=0 )
plot( d$year , d$temp , xlab="Year (CE)" , ylab="March temperature" , type="l" , lwd=2 , bty="l" , xlim=c(800,2015) , cex.axis=0.8 )
lines( d$year , d$temp_lower , col="gray" )
lines( d$year , d$temp_upper , col="gray" )
mtext( "March temperature reconstruction (95\% interval)" , adj=0 )

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

