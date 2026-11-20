# https://github.com/seantalts/stanli/blob/main/docs/from-cmdstanr.md

if ( FALSE ) {

library(rethinking)
data(chimpanzees)

# don't want any variables with NAs
# also recode condition to an index {1,0} -> {1,2}
d <- list( 
    pulled_left = chimpanzees$pulled_left ,
    prosoc_left = chimpanzees$prosoc_left ,
    condition = as.integer( 2 - chimpanzees$condition ) ,
    actor = as.integer( chimpanzees$actor ) ,
    blockid = as.integer( chimpanzees$block )
)

m3 <- ulam(
    alist(
        pulled_left ~ bernoulli(theta),
        logit(theta) <- a + aj[actor] + ak[blockid] + bp[condition]*prosoc_left,
        aj[actor] ~ normal( 0 , sigma_actor ),
        ak[blockid] ~ normal( 0 , sigma_block ),
        a ~ dnorm(0,4),
        bp[condition] ~ dnorm(0,1),
        sigma_actor ~ exponential(1),
        sigma_block ~ exponential(1)
    ) ,
    data=d, chains=2 , cores=2 , sample=TRUE )

m3s <- ulam(
    alist(
        pulled_left ~ bernoulli(theta),
        logit(theta) <- a + aj[actor] + ak[blockid] + bp[condition]*prosoc_left,
        aj[actor] ~ normal( 0 , sigma_actor ),
        ak[blockid] ~ normal( 0 , sigma_block ),
        a ~ dnorm(0,4),
        bp[condition] ~ dnorm(0,1),
        sigma_actor ~ exponential(1),
        sigma_block ~ exponential(1)
    ) ,
    data=d, chains=2 , cores=2 , sample=TRUE , stanli=TRUE )

mm <- attr(m3s,"stanlifit")
precis(m3s,3)
p <- extract.samples(m3s)


m4s <- ulam(
    alist(
        # likeliood
        pulled_left ~ bernoulli(theta),

        # linear models
        logit(theta) <- A + BP*prosoc_left,
        A <- a + v[actor,1],
        BP <- bp + v[actor,condition+1],

        # adaptive prior
        vector[3]: v[actor] ~ multi_normal( 0 , Rho_actor , sigma_actor ),

        # fixed priors
        c(a,bp) ~ normal(0,1),
        sigma_actor ~ exponential(1),
        Rho_actor ~ lkjcorr(4)
    ) , data=d , chains=3 , cores=3 , stanli=TRUE )

m5s <- ulam(
    alist(
        # likeliood
        pulled_left ~ bernoulli(theta),

        # linear models
        logit(theta) <- A + BP*prosoc_left,
        A <- a + v[actor,1],
        BP <- bp + v[actor,condition+1],

        # adaptive prior
        matrix[actor,3]: v <- compose_noncentered( sigma_actor , L_Rho_actor , z ),
        matrix[3,actor]: z ~ normal( 0 , 1 ),

        # fixed priors
        c(a,bp) ~ normal(0,1),
        vector[3]: sigma_actor ~ exponential(1),
        cholesky_factor_corr[3]: L_Rho_actor ~ lkj_corr_cholesky( 4 )
    ) , data=d , chains=3 , cores=3 , stanli=TRUE )

m6s <- ulam(
    alist(
        # likeliood
        pulled_left ~ bernoulli(theta),

        # linear models
        logit(theta) <- A + BP*prosoc_left,
        A <- a + v[actor,1],
        BP <- bp + v[actor,condition+1],

        # adaptive prior
        matrix[actor,3]: v <- t(diag_pre_multiply( sigma_actor , L_Rho_actor ) * z),
        matrix[3,actor]: z ~ normal( 0 , 1 ),

        # fixed priors
        c(a,bp,bpc) ~ normal(0,1),
        vector[3]: sigma_actor ~ exponential(1),
        cholesky_factor_corr[3]: L_Rho_actor ~ lkj_corr_cholesky( 4 )
    ) , data=d , chains=3 , cores=3 , stanli=TRUE )

m6 <- ulam(
    alist(
        # likeliood
        pulled_left ~ bernoulli(theta),

        # linear models
        logit(theta) <- A + BP*prosoc_left,
        A <- a + v[actor,1],
        BP <- bp + v[actor,condition+1],

        # adaptive prior
        matrix[actor,3]: v <- t(diag_pre_multiply( sigma_actor , L_Rho_actor ) * z),
        matrix[3,actor]: z ~ normal( 0 , 1 ),

        # fixed priors
        c(a,bp,bpc) ~ normal(0,1),
        vector[3]: sigma_actor ~ exponential(1),
        cholesky_factor_corr[3]: L_Rho_actor ~ lkj_corr_cholesky( 4 )
    ) , data=d , chains=3 , cores=3 , stanli=FALSE )

}
