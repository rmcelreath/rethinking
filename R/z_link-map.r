# function to compute value of linear models from map fit, over samples

#' Computes inverse-link linear model values for \code{map} and \code{map2stan}
#' samples.
#'
#' This function computes the value of each linear model at each sample for
#' each case in the data. Inverse link functions are applied, so that for
#' example a logit link linear model produces probabilities, using the logistic
#' transform.
#'
#' This function is used internally by \code{\link{WAIC}}, \code{\link{sim}},
#' \code{\link{postcheck}}, and \code{\link{ensemble}}.
#'
#' It is possible to replace components of the posterior distribution with
#' simulated values. The \code{replace} argument should be a named \code{list}
#' with replacement values. This is useful for marginalizing over varying
#' effects. See the examples below for an example in which varying intercepts
#' are marginalized this way.
#'
#' @name link-methods
#' @aliases link link-methods link,map2stan-method link,map-method
#' link,lm-method
#' @docType methods
#' @param fit Object of class \code{map} or \code{map2stan}
#' @param data Optional list of data to compute predictions over. When missing,
#' uses data found inside fit object.
#' @param n Number of samples to use
#' @param post Optional samples from posterior. When missing, \code{link}
#' extracts the samples using \code{\link{extract.samples}}.
#' @param refresh Refresh interval for progress display. Set to
#' \code{refresh=0} to suppress display.
#' @param replace Optional named list of samples to replace inside posterior
#' samples. See examples.
#' @param flatten When \code{TRUE}, removes linear model names from result
#' @param ... Other parameters to pass to someone
#' @author Richard McElreath
#' @seealso \code{\link{map}}, \code{\link{map2stan}}, \code{\link{sim}},
#' \code{\link{ensemble}}, \code{\link{postcheck}}
#' @export
#' @examples
#'
#' \dontrun{
#' library(rethinking)
#' data(chimpanzees)
#' d <- chimpanzees
#' d$recipient <- NULL     # get rid of NAs
#'
#' # model 4 from chapter 12 of the book
#' m12.4 <- map2stan(
#'     alist(
#'         pulled_left ~ dbinom( 1 , p ) ,
#'         logit(p) <- a + a_actor[actor] + (bp + bpC*condition)*prosoc_left ,
#'         a_actor[actor] ~ dnorm( 0 , sigma_actor ),
#'         a ~ dnorm(0,10),
#'         bp ~ dnorm(0,10),
#'         bpC ~ dnorm(0,10),
#'         sigma_actor ~ dcauchy(0,1)
#'     ) ,
#'     data=d , warmup=1000 , iter=4000 , chains=4 )
#'
#' # posterior predictions for a particular actor
#' chimp <- 2
#' d.pred <- list(
#'     prosoc_left = c(0,1,0,1),   # right/left/right/left
#'     condition = c(0,0,1,1),     # control/control/partner/partner
#'     actor = rep(chimp,4)
#' )
#' link.m12.4 <- link( m12.4 , data=d.pred )
#' apply( link.m12.4 , 2 , mean )
#' apply( link.m12.4 , 2 , PI )
#'
#' # posterior predictions marginal of actor
#' # here we replace the varying intercepts samples
#' #   with simulated values
#'
#' # replace varying intercept samples with simulations
#' post <- extract.samples(m12.4)
#' a_actor_sims <- rnorm(7000,0,post$sigma_actor)
#' a_actor_sims <- matrix(a_actor_sims,1000,7)
#'
#' # fire up link
#' # note use of replace list
#' link.m12.4 <- link( m12.4 , n=1000 , data=d.pred ,
#'     replace=list(a_actor=a_actor_sims) )
#'
#' # summarize
#' apply( link.m12.4 , 2 , mean )
#' apply( link.m12.4 , 2 , PI )
#' }
#'
#' @export
setGeneric("link",
function( fit , data , n=1000 , ... ) {
    print(class(fit))
}
)

#' @export
setMethod("link", "map",
function( fit , data , n=1000 , post , refresh=0.1 , replace=list() , flatten=TRUE , ... ) {

    if ( ! inherits(fit, "map")) stop("Requires map fit")
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
