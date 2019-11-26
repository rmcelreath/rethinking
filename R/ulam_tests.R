#######################################
#######################################
# ulam test code
if ( FALSE ) {

    library(rethinking)

    # compare hash of output to target hash
    check_hash <- function( x , h ) {
        require(openssl)
        if ( md5(x) != h )
            stop( concat( "Hashes not identical:" , paste( match.call() , collapse=" , " ) ) )
    }

    # student t test

    library(rethinking)
    data(Howell1)
    d <- Howell1
    d2 <- d[ d$age >= 18 , ]

    z <- ulam(
      alist(
        height ~ student_t( nu, mu , sigma ) ,
        mu <- a + b*weight ,
        a ~ normal( 178 , 20 ) ,
        b ~ normal( 0 , 1 ) ,
        sigma ~ cauchy( 0 , 2),
        nu ~ gamma(2, 0.1)
      ) , warmup=1000,iter=2000,chains=2,cores=2,
    data=d2 )

    s <- sim(z)

    # dbeta2 test

    y <- rbeta2( 100 , 0.3 , 3 )

    z <- ulam(
        alist(
            y ~ dbeta2( p , theta ),
            logit(p) <- a,
            a ~ normal(0,1.5),
            theta ~ exponential(1)
        ), data=list(y=y) , sample=TRUE )

    # test "no prior" error
    z <- ulam(
        alist(
            y ~ dbeta2( p , theta ),
            logit(p) <- a,
            #a ~ normal(0,1.5),
            theta ~ exponential(1)
        ), data=list(y=y) , sample=TRUE )

    z <- ulam(
        alist(
            y ~ dbeta2( p , theta ),
            #logit(p) <- a,
            #a ~ normal(0,1.5),
            theta ~ exponential(1)
        ), data=list(y=y) , sample=TRUE )

    # binom tests

    library(rethinking)
    data( UCBadmit )
    UCBadmit$male <- as.integer( ifelse( UCBadmit$applicant.gender=="male" , 1 , 0 ) )
    UCBadmit$dept <- rep( 1:6 , each=2 )
    UCBadmit$ag <- UCBadmit$applicant.gender
    UCBadmit$applicant.gender <- NULL

    # quap to ulam
    z <- quap(
        alist(
            admit ~ dbinom(applications,p),
            logit(p) <- a,
            a ~ dnorm(0,4)
        ),
        data=UCBadmit )
    
    zz <- ulam( z , sample=FALSE )

    # binomial
    z <- ulam(
        alist(
            admit ~ binomial(applications,p),
            logit(p) <- a,
            a ~ normal(0,4)
        ),
        data=UCBadmit , sample=FALSE , log_lik=TRUE )
    check_hash( stancode(z) , "2b8d92e0c80cffd990a526556c72308e" )

    # abitrary constraints
    z <- ulam(
        alist(
            admit ~ binomial(applications,p),
            logit(p) <- a,
            a ~ normal(0,4)
        ),
        constraints=list(a="lower=0"),
        data=UCBadmit , sample=TRUE , log_lik=FALSE , chains=4 )

    z <- ulam(
        alist(
            admit ~ binomial(applications,p),
            logit(p) <- a,
            a ~ normal(0,4)
        ),
        data=UCBadmit , sample=FALSE , log_lik=TRUE )

    # test to make sure constant '1' in binomial here isn't [i]'d in log_lik code
    data( chimpanzees )
    z <- ulam(
        alist(
            y ~ binomial(1,p),
            logit(p) <- a,
            a ~ normal(0,4)
        ),
        data=list(y = chimpanzees$pulled_left) , sample=FALSE , log_lik=TRUE )

    z2 <- ulam(
        alist(
            admit ~ binomial(applications,p),
            logit(p) <- a + b*male,
            a ~ normal(0,4),
            b ~ normal(0,1)
        ),
        data=UCBadmit , log_lik=TRUE , sample=FALSE )

    check_hash( stancode(z2) , "30e2efaf60da86c26c12d1546be851c2" )

    # continuous missing data

    # automated merging and NA indexing
    # merge_missing
    UCBadmit$x <- rnorm(12)
    UCBadmit$x[1:2] <- NA
    UCBadmit$x[2] <- NA

    z2b2 <- ulam(
        alist(
            admit ~ binomial(applications,p),
            logit(p) <- a + b*male + bx*x_merge,
            x_merge ~ normal( 0 , 1 ),
            x_merge <- merge_missing( x , x_impute ),
            a ~ normal(0,4),
            b ~ normal(0,1),
            bx ~ normal(0,1)
        ),
        data=UCBadmit , sample=TRUE , log_lik=FALSE )
    check_hash( stancode(z2b) , "9340600d32939e0ab63999a8aefc7663" )

    # automated version
    # needs to build the merge code
    z2b_auto <- ulam(
        alist(
            admit ~ binomial(applications,p),
            logit(p) <- a + b*male + bx*x,
            x ~ normal( 0 , 1 ),
            a ~ normal(0,4),
            b ~ normal(0,1),
            bx ~ normal(0,1)
        ),
        data=UCBadmit , sample=TRUE , log_lik=TRUE )
    check_hash( stancode(z2b_auto) , "9340600d32939e0ab63999a8aefc7663" )

    UCBadmit$x2 <- rnorm(12)
    UCBadmit$x2[5:7] <- NA

    z2c <- ulam(
        alist(
            admit ~ binomial(applications,p),
            logit(p) <- a + b*male + bx*x_merge + bx2*x2_merge,
            x_merge ~ normal( 0 , 1 ),
            x2_merge ~ normal( 0 , 1 ),
            x_merge <- merge_missing( x , x_impute ),
            x2_merge <- merge_missing( x2 , x2_impute ),
            a ~ normal(0,4),
            b ~ normal(0,1),
            c(bx,bx2) ~ normal(0,1)
        ),
        data=UCBadmit , sample=FALSE )

    # discrete? need to mix over different terms
    library(rethinking)
    data( UCBadmit )
    UCBadmit$male <- as.integer( ifelse( UCBadmit$applicant.gender=="male" , 1 , 0 ) )
    UCBadmit$dept <- rep( 1:6 , each=2 )
    UCBadmit$applicant.gender <- NULL
    UCBadmit$male2 <- UCBadmit$male
    UCBadmit$male2[1:2] <- (-1) # missingness code
    UCBadmit$male2 <- as.integer(UCBadmit$male2)

    z <- ulam(
        alist(
            admit|male2==-1 ~ mixture(  
                phi_male , 
                binomial( applications , p_m1 ) , 
                binomial( applications , p_m0) ),
            admit|male2>-1 ~ binomial( applications , p ),
            logit(p) <- a[dept] + b*male2,
            logit(p_m1) <- a[dept] + b*1,
            logit(p_m0) <- a[dept] + b*0,
            male2|male2>-1 ~ bernoulli( phi_male ),
            phi_male ~ beta(2,2),
            a[dept] ~ normal(0,4),
            b ~ normal(0,1)
        ),
        data=UCBadmit , sample=FALSE , log_lik=TRUE )

    z <- ulam(
        alist(
            admit|male2==-1 ~ mixture(  
                phi_male , 
                binomial( applications , p_m2 ) ,
                binomial( applications , p_m1 ) , 
                binomial( applications , p_m0) ),
            admit|male2>-1 ~ binomial( applications , p ),
            logit(p) <- a[dept] + b*male2,
            logit(p_m1) <- a[dept] + b*1,
            logit(p_m0) <- a[dept] + b*0,
            male2|male2>-1 ~ bernoulli( phi_male ),
            phi_male ~ beta(2,2),
            a[dept] ~ normal(0,4),
            b ~ normal(0,1)
        ),
        data=UCBadmit , sample=FALSE , log_lik=TRUE )

    z2d <- ulam(
        alist(
            admit|male2==-1 ~ custom( log_mix( 
                phi_male , 
                binomial_lpmf(admit|applications,p_m1) , 
                binomial_lpmf(admit|applications,p_m0) ) ),
            admit|male2>-1 ~ binomial( applications , p ),
            logit(p) <- a[dept] + b*male2,
            logit(p_m1) <- a[dept] + b*1,
            logit(p_m0) <- a[dept] + b*0,
            male2|male2>-1 ~ bernoulli( phi_male ),
            phi_male ~ beta(2,2),
            a[dept] ~ normal(0,4),
            b ~ normal(0,1)
        ),
        data=UCBadmit , sample=FALSE , log_lik=TRUE )

    check_hash( stancode(z2d) , "e4ac93844f2c0e346fde546dbe7e8906" )

    # log_sum_exp form
    z2e <- ulam(
        alist(
            admit|male2==-1 ~ custom( 
                log_sum_exp( 
                    log(phi_male) + binomial_lpmf(admit|applications,p_m1) , 
                    log1m(phi_male) + binomial_lpmf(admit|applications,p_m0) ) ),
            admit|male2>-1 ~ binomial( applications , p ),
            logit(p) <- a + b*male2,
            logit(p_m1) <- a + b*1,
            logit(p_m0) <- a + b*0,
            male2|male2>-1 ~ bernoulli( phi_male ),
            phi_male ~ beta(2,2),
            a ~ normal(0,4),
            b ~ normal(0,1)
        ),
        data=UCBadmit , sample=FALSE )

    # bracketed intercepts
    z3 <- ulam(
        alist(
            admit ~ binomial(applications,p),
            logit(p) <- a[dept] + b*male,
            a[dept] ~ normal(0,4),
            b ~ normal(0,1)
        ),
        data=UCBadmit , sample=TRUE , log_lik=TRUE )

    z3r <- ulam(
        alist(
            admit ~ dbinom(applications,p),
            logit(p) <- a[dept] + b*male,
            a[dept] ~ dnorm(0,4),
            b ~ dnorm(0,1)
        ),
        data=UCBadmit , sample=FALSE )

    # gaussian outcome
    UCBadmit$y <- UCBadmit$admit/UCBadmit$applications
    z3b <- ulam(
        alist(
            y ~ normal(p,sigma),
            p <- a[dept] + b*male,
            a[dept] ~ normal(0,4),
            b ~ normal(0,1),
            sigma ~ exponential(1)
        ),
        data=UCBadmit )

    z4 <- ulam(
        alist(
            admit ~ binomial(applications,p),
            logit(p) <- a[dept] + b*male,
            a[dept] ~ normal( abar , sigma ),
            abar ~ normal( 0 , 4 ),
            sigma ~ normal(0,1),
            b ~ normal(0,1)
        ),
        data=UCBadmit )

    z5 <- ulam(
        alist(
            admit ~ binomial(applications,p),
            logit(p) <- a[dept] + b[dept]*male,
            c( a , b )[dept] ~ multi_normal( c(abar,bbar) , Rho , sigma ),
            abar ~ normal( 0 , 4 ),
            bbar ~ normal(0,1),
            sigma ~ half_normal(0,1),
            Rho ~ lkjcorr(2)
        ),
        data=UCBadmit , sample=FALSE , log_lik=TRUE )

    z6 <- ulam(
        alist(
            admit ~ binomial(applications,p),
            logit(p) <- v[dept,1] + v[dept,2]*male,
            vector[2]:v[dept] ~ multi_normal( c(abar,bbar) , Rho , sigma ),
            abar ~ normal( 0 , 4 ),
            bbar ~ normal(0,1),
            sigma ~ half_normal(0,1),
            Rho ~ lkjcorr(2)
        ),
        data=UCBadmit , sample=TRUE , log_lik=TRUE )

    z6b <- ulam(
        alist(
            admit ~ binomial(applications,p),
            logit(p) <- abar + v[dept,1] + ( bbar + v[dept,2] )*male,
            vector[2]:v[dept] ~ multi_normal( 0 , Rho , sigma ),
            abar ~ normal( 0 , 4 ),
            bbar ~ normal(0,1),
            sigma ~ half_normal(0,1),
            Rho ~ lkjcorr(2)
        ),
        data=UCBadmit , sample=TRUE )

    # version with second linear model for slope
    z6c <- ulam(
        alist(
            admit ~ binomial(applications,p),
            logit(p) <- abar + v[dept,1] + B*male,
            B <- bbar + v[dept,2],
            vector[2]:v[dept] ~ multi_normal( 0 , Rho , sigma ),
            abar ~ normal( 0 , 4 ),
            bbar ~ normal(0,1),
            sigma ~ half_normal(0,1),
            Rho ~ lkjcorr(2)
        ),
        data=UCBadmit )

    # version with TWO linear models for intercept, slope
    z6d <- ulam(
        alist(
            admit ~ binomial(applications,p),
            logit(p) <- A + B*male,
            A <- abar + v[dept,1],
            B <- bbar + v[dept,2],
            vector[2]:v[dept] ~ multi_normal( 0 , Rho , sigma ),
            abar ~ normal( 0 , 4 ),
            bbar ~ normal(0,1),
            sigma ~ half_normal(0,1),
            Rho ~ lkjcorr(2)
        ),
        data=UCBadmit , sample=TRUE )

    z7 <- ulam(
        alist(
            admit ~ binomial(applications,p),
            logit(p) <- v[dept,1] + v[dept,2]*male,
            vector[2]:v[dept] ~ multi_normal( c(abar,bbar) , Rho , sigma ),
            abar ~ normal( 0 , 4 ),
            bbar ~ normal(0,1),
            sigma[1] ~ half_normal(0,1),
            sigma[2] ~ half_normal(0,2),
            Rho ~ lkjcorr(2)
        ),
        data=UCBadmit )

    z7b <- ulam(
        alist(
            admit ~ binomial(applications,p),
            logit(p) <- v[dept,1] + v[dept,2]*male,
            vector[2]:v[dept] ~ multi_normal( v_mu , Rho , sigma ),
            vector[2]:v_mu ~ normal(0,1),
            sigma[1] ~ half_normal(0,1),
            sigma[2] ~ half_normal(0,2),
            Rho ~ lkjcorr(2)
        ),
        data=UCBadmit )

    z7c <- ulam(
        alist(
            admit ~ binomial(applications,p),
            logit(p) <- v[dept,1] + v[dept,2]*male,
            vector[2]:v[dept] ~ multi_normal( v_mu , Rho , sigma ),
            vector[2]:v_mu[[1]] ~ normal(0,4),
            vector[2]:v_mu[[2]] ~ normal(0,1),
            sigma[1] ~ half_normal(0,1),
            sigma[2] ~ half_normal(0,2),
            Rho ~ lkjcorr(2)
        ),
        data=UCBadmit )

    # cholesky form

    # multi_normal_NC needs to edit formula with macro
    # DOES NOT YET WORK
    z8p <- ulam(
        alist(
            admit ~ binomial(applications,p),
            logit(p) <- v_mu[1] + v[dept,1] + (v_mu[2] + v[dept,2])*male,
            vector[2]:v[dept] ~ multi_normal_NC( Rho , sigma ),
            vector[2]: v_mu[[1]] ~ normal(0,4),
            vector[2]: v_mu[[2]] ~ normal(0,1),
            sigma ~ half_normal(0,1),
            Rho ~ lkjcorr( 2 )
        ),
        data=UCBadmit , chains=4 )

    z8 <- ulam(
        alist(
            admit ~ binomial(applications,p),
            logit(p) <- v_mu[1] + v[dept,1] + (v_mu[2] + v[dept,2])*male,
            matrix[dept,2]: v <- t(diag_pre_multiply( sigma , L_Rho ) * z),
            matrix[2,dept]: z ~ normal( 0 , 1 ),
            vector[2]: v_mu[[1]] ~ normal(0,4),
            vector[2]: v_mu[[2]] ~ normal(0,1),
            vector[2]: sigma ~ half_normal(0,1),
            cholesky_factor_corr[2]: L_Rho ~ lkj_corr_cholesky( 2 )
        ),
        data=UCBadmit , chains=4 )

    # version that saves the v vector by adding also to gen quants
    # note "save>" syntax
    z8c <- ulam(
        alist(
            admit ~ binomial(applications,p),
            logit(p) <- v_mu[1] + v[dept,1] + (v_mu[2] + v[dept,2])*male,
            save> matrix[dept,2]: v <- t(diag_pre_multiply( sigma , L_Rho ) * z),
            matrix[2,dept]: z ~ normal( 0 , 1 ),
            vector[2]: v_mu[[1]] ~ normal(0,4),
            vector[2]: v_mu[[2]] ~ normal(0,1),
            vector[2]: sigma ~ half_normal(0,1),
            cholesky_factor_corr[2]: L_Rho ~ lkj_corr_cholesky( 2 )
        ),
        data=UCBadmit )

    z8b <- ulam(
        alist(
            admit ~ binomial(applications,p),
            logit(p) <- v_mu[1] + v[dept,1] + (v_mu[2] + v[dept,2])*male,
            matrix[dept,2]: v <- compose_noncentered( sigma , L_Rho , z ),
            matrix[2,dept]: z ~ normal( 0 , 1 ),
            vector[2]: v_mu[[1]] ~ normal(0,4),
            vector[2]: v_mu[[2]] ~ normal(0,1),
            vector[2]: sigma ~ half_normal(0,1),
            cholesky_factor_corr[2]: L_Rho ~ lkj_corr_cholesky( 2 )
        ),
        data=UCBadmit )

    # multivariate model
    y <- rmvnorm( 100 , mean=c(1,2) , sigma=matrix(c(1,0.7,0.7,1),2,2) )
    z9 <- ulam(
        alist(
            y ~ multi_normal( MU , Rho , sigma ),
            vector[2]: MU ~ normal(0,1),
            corr_matrix[2]: Rho ~ lkj_corr(2),
            vector[2]: sigma ~ half_normal(0,1)
        ),
        data=list(y=y) , sample=FALSE )

    z9b <- ulam(
        alist(
            c(y1,y2) ~ multi_normal( MU , Rho , sigma ),
            vector[2]: MU ~ normal(0,1),
            corr_matrix[2]: Rho ~ lkj_corr(2),
            vector[2]: sigma ~ half_normal(0,1)
        ),
        data=list(y1=y[,1],y2=y[,2]) , sample=FALSE )

    # ordered logit test
    data(Trolley)
    dat_trolley <- list(
        response = Trolley$response,
        action = Trolley$action,
        intention = Trolley$intention,
        contact = Trolley$contact,
        id = coerce_index(Trolley$id)
    )

    z10 <- ulam(
        alist(
            response ~ ordered_logistic( theta , cutpoints ),
            theta <- ba*action + bi*intention,
            c(ba,bi) ~ normal(0,1),
            ordered[6]: cutpoints ~ normal(0,4)
        ),
        data=dat_trolley , sample=FALSE )

    z10b <- ulam(
        alist(
            response ~ ordered_logistic( theta , cutpoints ),
            theta <- ba*action + bi*intention,
            c(ba,bi) ~ normal(0,1),
            cutpoints ~ normal(0,4)
        ),
        data=dat_trolley , sample=FALSE )

    z10 <- ulam(
        alist(
            response ~ categorical_logit( theta ),
            theta <- ba*action + bi*intention,
            c(ba,bi) ~ normal(0,1)
        ),
        data=dat_trolley , sample=FALSE )

    # zero-inflated poisson
    # gen data first - example from text
    prob_drink <- 0.2 # 20% of days
    rate_work <- 1    # average 1 manuscript per day
    N <- 365
    drink <- rbinom( N , 1 , prob_drink )
    y <- as.integer( (1-drink)*rpois( N , rate_work ) )

    z11 <- ulam(
        alist(
            y ~ zipoisson( p , lambda ),
            logit(p) <- ap,
            log(lambda) <- al,
            ap ~ dnorm(0,1),
            al ~ dnorm(0,10)
        ) ,
        data=list(y=y) )

    z11b <- ulam(
        alist(
            y|y==0 ~ custom( log_mix( p , 0 , poisson_lpmf(0|lambda) ) ),
            y|y>0 ~ custom( log1m(p) + poisson_lpmf(y|lambda) ),
            logit(p) <- ap,
            log(lambda) <- al,
            ap ~ dnorm(0,1),
            al ~ dnorm(0,10)
        ) ,
        data=list(y=y) )

    x <- rnorm( N )
    z11c <- ulam(
        alist(
            y|y==0 ~ custom( log_mix( p , 0 , poisson_lpmf(0|lambda) ) ),
            y|y>0 ~ custom( log1m(p) + poisson_lpmf(y|lambda) ),
            logit(p) <- ap,
            log(lambda) <- al + bl*x,
            ap ~ dnorm(0,1),
            al ~ dnorm(0,10),
            bl ~ normal(0,1)
        ) ,
        data=list(y=y,x=x) , sample=FALSE )

    # GP?
    library(rethinking)
    data(Kline2)
    d <- Kline2
    data(islandsDistMatrix)
    d$society <- 1:10
    dat <- list(
        y=d$total_tools,
        society=d$society,
        log_pop = log(d$population),
        Dmat=islandsDistMatrix
    )

    z12 <- ulam(
        alist(
            y ~ poisson( mu ),
            log(mu) <- a + b*log_pop,
            a ~ normal(0,10),
            b ~ normal(0,1),
            etasq ~ exponential(1),
            rhosq ~ exponential(1)
        ),
        data=dat , sample=FALSE )

    z12b <- ulam(
        alist(
            y ~ poisson( mu ),
            log(mu) <- a + aj[society] + b*log_pop,
            a ~ normal(0,10),
            b ~ normal(0,1),
            vector[10]: aj ~ multi_normal( 0 , SIGMA ),
            matrix[10,10]: SIGMA <- cov_GPL2( Dmat , etasq , rhosq , 0.01 ),
            etasq ~ exponential(1),
            rhosq ~ exponential(1)
        ),
        data=dat , sample=TRUE , chains=4 )

    # Mark-recapture - Cormack-Jolly-Seber
    # individual capture history matrix
    # individuals in rows, surveys in columns

    # sim data
    K <- 5;                      # captures
    p <- runif(K, 0.7, 0.8);      # capture at k
    phi <- runif(K-1, 0.8, 0.9);  # survive k to k+1

    I <- 100;                    # individual
    X <- matrix(0,I,K);           # capture of individual
    for (i in 1:I) {
      alive <- 1;
      for (k in 1:K) {
        if (alive)
          X[i,k] <- rbinom(1,1,p[k]);
        if (alive && (k < K))
          alive <- rbinom(1,1,phi[k]);
      }
    }

    dat <- list( y = X )

    z13 <- ulam(
        alist(
            y ~ CJS_capture_history( p , phi ),
            vector[5]: p ~ beta( 2 , 2 ),
            vector[4]: phi ~ beta( 2 , 2 )
        ),
        data=dat , sample=FALSE )

}
