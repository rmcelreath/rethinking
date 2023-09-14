# test cmdstanr installation through ulam

test_ulam <- function() {
    m <- ulam(
        alist(
            x ~ normal(0,1),
            y ~ normal(1,1)
        ), data=list(a=1) )
    precis(m)
}
