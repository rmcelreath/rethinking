# sim -- simulates observations from fit map and map2stan models

setGeneric("sim",
function( fit , data , n=1000 , ... ) {
    predict(fit)
}
)

setMethod("sim", "map",
function( fit , data , n=1000 , post , ... ) {
    
    ########################################
    # check arguments
    if ( class(fit)!="map" ) stop("Requires map fit")
    if ( missing(data) ) data <- fit@data
    
    if ( missing(post) ) 
        post <- extract.samples(fit,n=n)
    else
        n <- length(post[[1]])
    
    # get linear model values from link
    # use our posterior samples, so later parameters have right correlation structure
    # don't flatten result, so we end up with a named list, even if only one element
    pred <- link( fit , data=data , n=n , post=post , flatten=FALSE , ... )
    
    # extract likelihood, assuming it is first element of formula
    lik <- flist_untag(fit@formula)[[1]]
    # discover outcome
    outcome <- as.character(lik[[2]])
    # discover likelihood function
    flik <- as.character(lik[[3]][[1]])
    # get simulation partner function
    rlik <- flik
    substr( rlik , 1 , 1 ) <- "r"
    # pull out parameters in likelihood
    pars <- vector(mode="list",length=length(lik[[3]])-1)
    for ( i in 1:length(pars) ) {
        pars[[i]] <- lik[[3]][[i+1]]
    }
    pars <- paste( pars , collapse=" , " )
    # build expression to evaluate
    n_cases <- length(data[[1]])
    xeval <- paste( rlik , "(" , n_cases , "," , pars , ")" , collapse="" )
    
    # simulate outcomes
    sim_out <- matrix(NA,nrow=n,ncol=n_cases)
    for ( s in 1:n ) {
        # build environment
        elm <- pred[[1]][s,] # extract s-th sample case calculations
        elm <- list(elm)
        names(elm)[1] <- names(pred)[1]
        e <- list( as.list(data) , as.list(post[s,]) , as.list(elm) )
        e <- unlist( e , recursive=FALSE )
        # evaluate
        sim_out[s,] <- eval(parse(text=xeval),envir=e)
    }
    
    return(sim_out)
}
)

# EXAMPLES/TEST
if ( FALSE ) {

library(rethinking)
data(cars)

fit <- map(
    alist(
        dist ~ dnorm(mu,sigma),
        mu <- a + b*speed
    ),
    data=cars,
    start=list(a=30,b=0,sigma=5)
)

pred <- link(fit,data=list(speed=0:30),flatten=FALSE)

sim.dist <- sim(fit,data=list(speed=0:30))

} #EXAMPLES
