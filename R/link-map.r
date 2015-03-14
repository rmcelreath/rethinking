# function to compute value of linear models from map fit, over samples

setGeneric("link",
function( fit , data , n=1000 , ... ) {
    print(class(fit))
}
)

setMethod("link", "map",
function( fit , data , n=1000 , post , refresh=0.1 , replace=list() , flatten=TRUE , ... ) {
    
    if ( class(fit)!="map" ) stop("Requires map fit")
    if ( missing(data) ) {
        data <- fit@data
    } else {
        # make sure all variables same length
        # weird vectorization errors otherwise
        #data <- as.data.frame(data)
    }
    
    if ( missing(post) ) 
        post <- extract.samples(fit,n=n)
    else {
        n <- dim(post[[1]])[1]
        if ( is.null(n) ) n <- length(post[[1]])
    }
    
    # replace with any elements of replace list
    if ( length( replace ) > 0 ) {
        for ( i in 1:length(replace) ) {
            post[[ names(replace)[i] ]] <- replace[[i]]
        }
    }
    
    nlm <- length(fit@links)
    f_do_lm <- TRUE
    if ( nlm==0 ) {
        # no linear models!
        # so need to flag to skip lm insertions into eval environment for ll
        nlm <- 1
        f_do_lm <- FALSE
    }
    
    link_out <- vector(mode="list",length=nlm)
    
    # for each linear model, compute value for each sample
    for ( i in 1:nlm ) {
        ref_inc <- floor(n*refresh)
        ref_next <- ref_inc
        
        if ( f_do_lm==TRUE ) {
            parout <- fit@links[[i]][[1]]
            lm <- fit@links[[i]][[2]]
        } else {
            parout <- "ll"
            lm <- "0"
            # hacky solution -- find density function and insert whatever expression in typical link spot
            flik <- as.character(fit@formula[[1]][[3]][[1]])
            # mu for Gaussian
            if ( flik=="dnorm" ) lm <- as.character( fit@formula[[1]][[3]][[2]] )
            # p for binomial -- assume in third spot, after size
            if ( flik=="dbinom" ) lm <- as.character( fit@formula[[1]][[3]][[3]] )
            # lambda for poisson
            if ( flik=="dpois" ) lm <- as.character( fit@formula[[1]][[3]][[2]] )
        }
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
            init <- list() # holds one row of samples across all params
            for ( j in 1:length(post) ) {
                par_name <- names(post)[ j ]
                dims <- dim( post[[par_name]] )
                # scalar
                if ( is.null(dims) ) init[[par_name]] <- post[[par_name]][s]
                # vector
                if ( length(dims)==2 ) init[[par_name]] <- post[[par_name]][s,]
                # matrix
                if ( length(dims)==3 ) init[[par_name]] <- post[[par_name]][s,,]
            }#j
            e <- list( as.list(data) , as.list(init) )
            e <- unlist( e , recursive=FALSE )
            value[s,] <- eval(parse(text=lm),envir=e)
        }#s
        link_out[[i]] <- value
        names(link_out)[i] <- parout
    }
    
    if ( refresh>0 ) cat("\n")
    
    if ( flatten==TRUE )
        if ( length(link_out)==1 ) link_out <- link_out[[1]]
    
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
        logit(p) <- a + b*prosoc.left,
        c(a,b) ~ dnorm(0,1)
    ),
    data=chimpanzees,
    start=list(a=0,b=0)
)

pred <- link(fit)

pred2 <- link(fit,data=list(prosoc.left=0:1),flatten=FALSE)

sim.pulls <- sim(fit,data=list(prosoc.left=0:1))

fit2 <- map2stan(
    alist(
        pulled.left ~ dbinom( 1 , p ),
        logit(p) <- a + b*prosoc.left,
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
