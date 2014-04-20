# function to compute value of linear models from map fit, over samples

setGeneric("link",
function( fit , data , n=1000 , ... ) {
    print(class(fit))
}
)

setMethod("link", "map",
function( fit , data , n=1000 , probs=NULL , refresh=0.1 , replace=list() , ... ) {

    if ( class(fit)!="map" ) stop("Requires map fit")
    if ( missing(data) ) data <- fit@data
    
    post <- extract.samples(fit,n=n)
    
    nlm <- length(fit@links)
    
    link_out <- vector(mode="list",length=nlm)
    
    # for each linear model, compute value for each sample
    for ( i in 1:nlm ) {
        ref_inc <- floor(n*refresh)
        ref_next <- ref_inc
        
        parout <- fit@links[[i]][[1]]
        lm <- fit@links[[i]][[2]]
        # empty matrix to hold samples-by-cases values of linear model
        value <- matrix(NA,nrow=n,ncol=length(data[[1]]))
        # for each sample
        for ( s in 1:n ) {
            # refresh progress display
            if ( refresh > 0 ) {
                if ( s == ref_next ) {
                    msg <- paste( "[" , s , "/" , n , "]\r" , collapse=" " )
                    cat(msg)
                    ref_next <- s + ref_inc
                    if ( ref_next > n ) ref_next <- n
                }
            }
        
            # make environment
            e <- list( as.list(data) , as.list(post[s,]) )
            e <- unlist( e , recursive=FALSE )
            value[s,] <- eval(parse(text=lm),envir=e)
        }
        link_out[[i]] <- value
        names(link_out)[i] <- parout
    }
    
    if ( refresh>0 ) cat("\n")
    return(link_out)
}
)

# TESTS
if (FALSE) {

library(rethinking)
data(chimpanzees)

fit <- map(
    alist(
        pulled.left ~ dbinom( 1 , p ),
        logit(p) ~ a + b*prosoc.left,
        c(a,b) ~ dnorm(0,1)
    ),
    data=chimpanzees,
    start=list(a=0,b=0)
)

pred <- link(fit)

pred2 <- link(fit,data=list(prosoc.left=0:1))

fit2 <- map2stan(
    alist(
        pulled.left ~ dbinom( 1 , p ),
        logit(p) ~ a + b*prosoc.left,
        c(a,b) ~ dnorm(0,1)
    ),
    data=list(
        pulled.left=chimpanzees$pulled.left,
        prosoc.left=chimpanzees$prosoc.left
    ),
    start=list(a=0,b=0)
)

preds <- link(fit2)

preds2 <- link(fit2,data=list(prosoc_left=0:1))

}
