# Prior/Likelihood/Posterior plots, plp

plp <- function( mapfit , ... ) {
    
    
    if ( !(class(mapfit) %in% c("map","map2stan")) ) {
        stop( "Requires model fit of class map or map2stan" )
    }
    is_m2s <- FALSE
    if ( class(mapfit)=="map2stan" ) is_m2s <- TRUE
    
    post <- extract.samples(mapfit)
    form <- mapfit@formula
    
}


#EXAMPLES/TESTS
if (FALSE) {

library(rethinking)
data(chimpanzees)
d <- chimpanzees

flist4 <- alist(
    pulled.left ~ dbinom( 1 , p ),
    logit(p) <- a + b*prosoc.left ,
    c(a,b) ~ dnorm(0,1)
)

fit <- map( flist4 , data=d , start=list(a=0,b=0) )



}
