
sim.train.test <- function( N=20 , k=3 , rho=c(0.15,-0.4) , b_sigma=100 , DIC=FALSE , WAIC=FALSE, devbar=FALSE , devbarout=FALSE ) {
    #require(MASS)
    n_dim <- 1+length(rho)
    if ( n_dim < k ) n_dim <- k
    Rho <- diag(n_dim)
    for ( i in 1:length(rho) ) {
        Rho[1,i+1] <- rho[i]
    }
    Rho[lower.tri(Rho)] <- Rho[upper.tri(Rho)]
    X.train <- mvrnorm( n=N , mu=rep(0,n_dim) , Sigma=Rho )
    X.test <- mvrnorm( n=N , mu=rep(0,n_dim) , Sigma=Rho )
    mm.train <- matrix(1,nrow=N,ncol=1)
    bnames <- "a"
    if ( k > 1 ) {
        mm.train <- cbind( mm.train , X.train[,2:k] )
        bnames <- c( "a" , paste("b",1:(k-1),sep="") )
        pnames <- paste("b",1:(k-1),sep="")
        P <- paste( "c(" , paste(pnames,collapse=",") , ")" , sep="" , collapse="" )
        Pf <- paste( P , "~ dnorm(0,",b_sigma,")" )
    }
    B <- paste( "c(" , paste(bnames,collapse=",") , ")" , sep="" , collapse="" )
    d <- list( y=X.train[,1] , mm=mm.train , Bvec=B )
    #m <- lm.fit( mm.train , X.train[,1] )
    flist <- list( y ~ dnorm( mu , 1 ) , mu ~ 0 + mm %*% eval(parse(text=Bvec)) )
    start.list <- list(a=0)
    if ( k>1 ) {
        flist[[3]] <- eval(parse(text=Pf))
        for ( i in 2:k ) {
            ptext <- paste( "b" , i-1 , sep="" , collapse="" )
            start.list[[i]] <- 0
            names(start.list)[i] <- ptext
        }
    }
    m <- map( flist , data=d , start=start.list )
    dev.train <- as.numeric( deviance(m) )
    mm.test <- matrix(1,nrow=N,ncol=1)
    if ( k > 1 ) mm.test <- cbind( mm.test , X.test[,2:k] )
    dev.test <- (-2)*sum( dnorm( X.test[,1] , mm.test %*% coef(m) , 1 , TRUE ) )
    # result
    result <- c( dev.train , dev.test )
    # DIC and WAIC
    if ( DIC==TRUE ) result <- c( result , as.numeric(DIC(m)) )
    if ( WAIC==TRUE ) {
        n_samples <- 1e3
        l <- link(m,n=n_samples,refresh=0)
        lppd <- 0
        pD <- 0
        for ( i in 1:N ) {
            ll <- dnorm( X.train[i,1] , l[,i] , 1 , TRUE )
            lppd <- lppd + log_sum_exp(ll) - log(n_samples)
            pD <- pD + var(ll)
        }
        result <- c( result , -2*(lppd-pD) )
    }
    if ( devbar==TRUE ) {
        post <- extract.samples( m , n=1e3 )
        dev <- sapply( 1:nrow(post) , function(i) 
            (-2)*sum( dnorm( X.train[,1] , mm.train %*% as.numeric(post[i,]) , 1 , TRUE ) )
        )
        result <- c( result , mean(dev) )
    }
    if ( devbarout==TRUE ) {
        post <- extract.samples( m , n=1e3 )
        dev <- sapply( 1:nrow(post) , function(i) 
            (-2)*sum( dnorm( X.test[,1] , mm.test %*% as.numeric(post[i,]) , 1 , TRUE ) )
        )
        result <- c( result , mean(dev) )
    }
    # return
    return(result)
}
