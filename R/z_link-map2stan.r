# posterior predictions for map2stan model fits
# defaults to using original data

# to do:
# (*) fix dzipois vectorization hack --- need to fix dzipois itself

# new link function that doesn't invoke Stan
setMethod("link", "map2stan",
function( fit , data , n=1000 , post , refresh=0.1 , replace=list() , flatten=TRUE , ... ) {

    if ( class(fit)!="map2stan" ) stop("Requires map2stan fit")
    if ( missing(data) ) {
        data <- fit@data
    } else {
        # make sure all variables same length
        # weird vectorization errors otherwise
        #data <- as.data.frame(data)
    }
    
    # get samples from Stan fit
    if ( missing(post) )
        post <- extract.samples(fit,n=n)
    
    # check n and truncate accordingly
    # if ( n==0 )
    
    # replace with any elements of replace list
    if ( length( replace ) > 0 ) {
        for ( i in 1:length(replace) ) {
            post[[ names(replace)[i] ]] <- replace[[i]]
        }
    }
    
    lm <- fit@formula_parsed$lm
    lik <- fit@formula_parsed$lik
    n_lm <- length(lm)
    f_do_lm <- TRUE
    n_lik <- length(lik)
    lm_names <- c()
    lm_lik <- c()
    if ( n_lm>0 ) {
        for ( i in 1:n_lm ) {
            lm_names <- c( lm_names , lm[[i]]$parameter )
            if ( n_lik > 0 ) {
                # find corresponding likelihood distribution
                for ( j in 1:n_lik ) {
                    if ( lik[[j]]$N_name == lm[[i]]$N_name ) {
                        lm_lik <- c( lm_lik , lik[[j]]$likelihood )
                    }
                }#j
            }
        }#i
    } else {
        stop( "There appear to be no linear models here" )
        # no linear models!
        # so need to flag to skip lm insertions into eval environment for ll
        n_lm <- 1
        f_do_lm <- FALSE
    }
    
    # number of samples
    n_samples <- dim( post[[1]] )[1]
    if ( n == 0 ) n <- n_samples # special flag for all samples in fit
    if ( n_samples < n ) n <- n_samples
    
    lm_out <- list()
    liks <- fit@formula_parsed$likelihood
    n_cases_list <- list()
    for ( i in 1:n_lm ) {
        # find out how many cases by finding corresponding likelihood and counting cases of outcome in data
        # can't just use old number of cases (from fit), because might have new user data
        n_cases <- 1
        K <- 1 # max level for ordered outcome
        N_name <- lm[[i]]$N_name
        outcome <- NA
        for ( j in 1:length(liks) ) {
            if ( liks[[j]]$N_name == N_name ) outcome <- liks[[j]]$outcome
        }
        if ( !is.na(outcome) ) {
            n_cases <- length( data[[ outcome ]] )
            if ( n_cases > 0 )
                K <- max( fit@data[[ outcome ]] ) # get max level from original data
            else {
                # outcome not found, so just fill with 1s
                n_cases <- length( data[[1]] ) # just hope is informative
                data[[outcome]] <- rep(1,n_cases)
            }
        }
        n_cases_list[[i]] <- n_cases
        # empty array with right dims
        lm_out[[ lm_names[i] ]] <- array( 0 , dim=c( n , n_cases ) )
        
        if ( lm_lik[[i]]=="ordered_logistic" ) {
            # need another dimension in result, because outcome is a vector
            lm_out[[ lm_names[i] ]] <- array( 0 , dim=c( n , n_cases , K ) )
        }
    }
    
    #init <- fit@start
    init <- list()
    
    ###################
    # loop over samples and compute each case for each linear model
    
    # initialize refresh counter
    ref_inc <- floor(n*refresh)
    ref_next <- ref_inc
    
    # prep reused objects
    rhs <- list()
    for ( k in 1:n_lm ) {
        # ready linear model code
        if ( f_do_lm==TRUE ) {
            rhs0 <- fit@formula_parsed$lm[[k]]$RHS
            rhs0 <- gsub( "[i]" , "" , rhs0 , fixed=TRUE )
        } else {
            # no linear models
            # so use dummy linear model with appropriate likelihood parameter
            parout <- "ll"
            rhs0 <- "0"
            # hacky solution -- find density function and insert whatever expression in typical link spot
            flik <- as.character(liks[[1]][[3]][[1]])
            # mu for Gaussian
            if ( flik=="dnorm" ) rhs0 <- as.character(liks[[1]][[3]][[2]])
            # p for binomial -- assume in third spot, after size
            if ( flik=="dbinom" ) rhs0 <- as.character(liks[[1]][[3]][[3]])
            # lambda for poisson
            if ( flik=="dpois" ) rhs0 <- as.character(liks[[1]][[3]][[2]])
        }
        # check for imputed predictors
        if ( length(fit@formula_parsed$impute_bank) > 0 ) {
            # replace each imputed variable name with merged name
            for ( kk in 1:length(fit@formula_parsed$impute_bank) ) {
                name_original <- names(fit@formula_parsed$impute_bank)[kk]
                name_merge <- concat( name_original , "_merge" )
                # _merge name should already be in linear model
                #rhs0 <- gsub( name_merge , concat(name_merge,"[i,]") , rhs0 , fixed=TRUE )
                #rhs0 <- gsub( name_original , name_merge , rhs0 , fixed=TRUE )
                # construct merged matrix in posterior
                #   this "parameter" is constant in columns for observed
                #   but has samples in columns for imputed
                n_cases <- length( data[[ name_original ]] )
                var_merged <- matrix( NA , nrow=n , ncol=n_cases )
                name_impute <- concat( name_original , "_impute" )
                missingness <- fit@formula_parsed$impute_bank[[kk]]$missingness
                for ( a_case in 1:n_cases ) {
                    if ( a_case %in% missingness ) {
                        # insert column of samples
                        idx <- which( missingness==a_case )
                        if ( length(missingness)>1 )
                            var_merged[ , a_case ] <- post[[name_impute]][, idx ]
                        else
                            var_merged[ , a_case ] <- post[[name_impute]]
                    } else {
                        # insert column of observed values
                        var_merged[ , a_case ] <- rep( data[[ name_original ]][a_case] , n )
                    }
                }#a_case
                # insert matrix into posterior
                post[[name_merge]] <- var_merged
                # insert template into init list
                init[[name_merge]] <- rep(0,n_cases)
            }#kk
        }
        # store
        rhs[[k]] <- rhs0
    }#k
    
    # loop over samples
    for ( i in 1:n ) {
    
        # refresh progress display
        if ( refresh > 0 ) {
            if ( i == ref_next ) {
                msg <- paste( "[" , i , "/" , n , "]\r" , collapse=" " )
                cat(msg)
                ref_next <- i + ref_inc
                if ( ref_next > n ) ref_next <- n
            }
        }
        
        # build inits
        for ( j in 1:length(post) ) {
            par_name <- names(post)[ j ]
            dims <- dim( post[[par_name]] )
            # scalar
            if ( length(dims)==1 ) init[[par_name]] <- post[[par_name]][i]
            # vector
            if ( length(dims)==2 ) init[[par_name]] <- post[[par_name]][i,]
            # matrix
            if ( length(dims)==3 ) init[[par_name]] <- post[[par_name]][i,,]
        }#j
        
        # loop over linear models and compute by pushing samples through stanfit
        # pass through linear models in reverse order, so lower models can be available for higher ones
        for ( k in n_lm:1 ) {
        
            # ready environment
            e <- list( as.list(data) , as.list(init) )
            e <- unlist( e , recursive=FALSE )

            # evaluate
            r <- eval(parse(text=rhs[[k]]),envir=e)
            flink <- fit@formula_parsed$lm[[k]]$link
            if ( flink=="log" ) r <- exp(r)
            if ( flink=="logit" ) r <- logistic(r)
            
            # special processing for ordered logit
            # this needs to be sped up
            if ( !is.null(lm_lik) ) {
                if ( lm_lik[k] == "ordered_logistic" ) {
                    # need vector of probabilities as result
                    cuts <- init[['cutpoints']]
                    v <- predict_ordlogit( r , cuts )
                    #r <- t(v)
                    r <- v
                }
            }
            
            # store
            if ( lm_lik[k] == "ordered_logistic" ) {
                lm_out[[ lm_names[k] ]][i,,] <- r
            } else {
                lm_out[[ lm_names[k] ]][i,] <- r
            }
            
            # make linear models values available for next linear models
            init[[ lm_names[k] ]] <- r
            
        }#k
        
    }#i
    
    if ( refresh>0 ) cat("\n")
    
    if ( flatten==TRUE )
        if ( length(lm_out)==1 ) lm_out <- lm_out[[1]]
    
    return( lm_out )
}
)

predict_ordlogit <- function( phi , a ) {
    K <- 1:(length(a)+1)
    a <- c( as.numeric(a) , Inf )
    p <- sapply( K , function(k) logistic(a[k]-phi) )
    na <- c(-Inf,a)
    np <- sapply( K , function(k) logistic(na[k]-phi) )
    p <- p - np
    return(p)
}
# predict.ordlogit( c(-1,0,1) , seq(-2,2,length.out=6) )




# test code
if ( FALSE ) {

library(rethinking)
data(willowtitn)
d <- willowtitn

m <- map2stan(
    alist(
        n ~ dzipois( p , lambda ),
        log(lambda) <- an + bn*percforest ,
        logit(p) ~ ap + bp*percforest ,
        c(an,ap,bn,bp) ~ dnorm(0,10)
    ),
    data=d,
    start=list(ap=0,an=1,bn=0,bp=0)
)

pred <- link( m )

d.new <- list(
    n = 1:4,
    percforest = c(10,20,30,40),
    N = 4
)
pred <- link( m , data=d.new )

######

library(rethinking)
data(UCBadmit)
d <- UCBadmit
d$male <- ifelse( d$applicant.gender=="male" , 1 , 0 )
d$dept <- as.integer(d$dept)

m <- map2stan(
    alist(
        admit ~ dbinom( applications , p ),
        logit(p) <- a + ak + b*male,
        ak[dept] ~ dnorm(0,sigma),
        a ~ dnorm(0,10),
        b ~ dnorm(0,10),
        sigma ~ dcauchy(0,2.5)
    ) ,
    data=list(admit=d$admit,applications=d$applications,male=d$male,dept=d$dept) ,
    start=list(a=0,b=0,ak=rep(0,6),sigma=1)
)

pred <- link( m )

y.mean <- apply( pred$p , 2 , mean )
y.ci <- apply( pred$p , 2 , PCI )

plot( d$admit/d$applications , ylim=c(0,1) )
lines( 1:12 , y.mean )
shade( y.ci , 1:12 )

# now with random slopes
m2 <- map2stan(
    list(
        admit ~ dbinom( applications , p ),
        logit(p) ~ a + ak + (b+bk)*male,
        c(ak,bk)|dept ~ dmvnorm(0,Sigma),
        a ~ dnorm(0,10),
        b ~ dnorm(0,10),
        Sigma ~ inv_wishart(3,diag(2))
    ) ,
    data=list(admit=d$admit,applications=d$applications,male=d$male,dept=d$dept) ,
    start=list(a=0,b=0,ak=rep(0,6),bk=rep(0,6),Sigma=diag(2))
)

pred <- link.map2stan( m2 )

y.mean <- apply( pred$p , 2 , mean )
y.ci <- apply( pred$p , 2 , PCI )

plot( d$admit/d$applications , ylim=c(0,1) )
lines( 1:12 , y.mean )
shade( y.ci , 1:12 )


#####

library(rethinking)
data(cars)

m <- map2stan(
    list(
        y ~ dnorm( mu , sigma ),
        mu ~ a + b*x,
        a ~ dnorm(0,10),
        b ~ dnorm(0,10),
        sigma ~ dcauchy(0,1)
    ) ,
    data=list(y=cars$dist,x=cars$speed) ,
    start=list(a=mean(cars$dist),b=0,sigma=10)
)

mu <- link.map2stan( m )

y.mean <- apply( mu$mu , 2 , mean )
y.ci <- apply( mu$mu , 2 , PCI )

plot( cars$dist )
lines( 1:50 , y.mean )
shade( y.ci , 1:50 )

# new data

dn <- list( y=1:50 , x=seq(0,30,length.out=50) , N=50 )

mu <- link.map2stan( m , data=dn )

y.mean <- apply( mu$mu , 2 , mean )
y.ci <- apply( mu$mu , 2 , PCI )

plot( cars$speed , cars$dist )
lines( dn$x , y.mean )
for ( p in seq(0.1,0.99,length.out=6) ) {
    y.ci <- apply( mu$mu , 2 , PCI , prob=p )
    shade( y.ci , dn$x , col=col.alpha("slateblue",0.15) )
}

}

