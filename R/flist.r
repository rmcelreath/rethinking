# alist replacement that allows = in addition to <- in formulas
# doesn't work with link functions...not sure if this can be fixed

flist <- function(...) {
    L <- match.call() # named elements in language tree
    nn <- names(L)
    print(deparse(L))
    for ( i in 2:length(L) ) {
        if ( nn[i]!="" ) {
            # named argument, so convert = to <- and reformat
            the_name <- names(L)[i]
            the_arg <- deparse(L[[i]])
            new_arg <- concat( the_name , " <- " , the_arg )
            L[[i]] <- str2expression(new_arg)[[1]]
            names(L)[i] <- ""
        }
    }#i
    # replace flist with alist
    L[[1]] <- str2expression("alist")[[1]]
    return(eval(L))
}

if (FALSE) {

# test

x <- flist(
    y ~ normal(mu,sigma),
    mu = a + b*x
)

N <- 100
Y <- rnorm(N)
X <- rnorm(N)

m00 <- ulam(
    flist(
        Y ~ dnorm(mu,1),
        mu = a + b*X,
        a ~ dnorm(0,1),
        b ~ dnorm(0,1)
    ) , data=list(Y=Y,X=X) )

}