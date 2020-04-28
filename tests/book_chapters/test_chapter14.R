# test_dir('rethinking/tests/book_chapters',reporter="summary")
# test_dir('rethinking/tests/book_chapters',filter="chapter14")

context('chapter 14')
library(rethinking)

expect_equiv_eps <- function( x , y , eps=0.01 ) {
    expect_equivalent( x , y , tolerance=eps )
}


## R code 14.1
a <- 3.5            # average morning wait time
b <- (-1)           # average difference afternoon wait time
sigma_a <- 1        # std dev in intercepts
sigma_b <- 0.5      # std dev in slopes
rho <- (-0.7)       # correlation between intercepts and slopes

## R code 14.2
Mu <- c( a , b )

## R code 14.3
cov_ab <- sigma_a*sigma_b*rho
Sigma <- matrix( c(sigma_a^2,cov_ab,cov_ab,sigma_b^2) , ncol=2 )

## R code 14.6
N_cafes <- 20

## R code 14.7
library(MASS)
set.seed(5) # used to replicate example
vary_effects <- mvrnorm( N_cafes , Mu , Sigma )

## R code 14.8
a_cafe <- vary_effects[,1]
b_cafe <- vary_effects[,2]

## R code 14.10
set.seed(22)
N_visits <- 10
afternoon <- rep(0:1,N_visits*N_cafes/2)
cafe_id <- rep( 1:N_cafes , each=N_visits )
mu <- a_cafe[cafe_id] + b_cafe[cafe_id]*afternoon
sigma <- 0.5  # std dev within cafes
wait <- rnorm( N_visits*N_cafes , mu , sigma )
d <- data.frame( cafe=cafe_id , afternoon=afternoon , wait=wait )

## R code 14.12
test_that("R code 14.12",{
    set.seed(867530)
    expect_warning(
    m14.1 <- ulam(
        alist(
            wait ~ normal( mu , sigma ),
            mu <- a_cafe[cafe] + b_cafe[cafe]*afternoon,
            c(a_cafe,b_cafe)[cafe] ~ multi_normal( c(a,b) , Rho , sigma_cafe ),
            a ~ normal(5,2),
            b ~ normal(-1,0.5),
            sigma_cafe ~ exponential(1),
            sigma ~ exponential(1),
            Rho ~ lkj_corr(2)
        ) , data=d , chains=4 , cores=4 )
    )
    expect_equivalent( dim(precis(m14.1,3)) , c(49,6) )
})

## R code 14.18
library(rethinking)
data(chimpanzees)
d <- chimpanzees
d$block_id <- d$block
d$treatment <- 1L + d$prosoc_left + 2L*d$condition

dat <- list(
    L = d$pulled_left,
    tid = d$treatment,
    actor = d$actor,
    block_id = as.integer(d$block_id) )

test_that("R code 14.18",{
    set.seed(4387510)
    expect_warning( m14.2 <- ulam(
        alist(
            L ~ dbinom(1,p),
            logit(p) <- g[tid] + alpha[actor,tid] + beta[block_id,tid],

            # adaptive priors
            vector[4]:alpha[actor] ~ multi_normal(0,Rho_actor,sigma_actor),
            vector[4]:beta[block_id] ~ multi_normal(0,Rho_block,sigma_block),

            # fixed priors
            g[tid] ~ dnorm(0,1),
            sigma_actor ~ dexp(1),
            Rho_actor ~ dlkjcorr(4),
            sigma_block ~ dexp(1),
            Rho_block ~ dlkjcorr(4)
        ) , data=dat , chains=4 , cores=4 ) )
    expect_equivalent( dim(precis(m14.2,3)) , c(96,6) )
})

## R code 14.19
test_that("R code 14.19",{
    set.seed(4387510)
    expect_warning(
    m14.3 <- ulam(
        alist(
            L ~ binomial(1,p),
            logit(p) <- g[tid] + alpha[actor,tid] + beta[block_id,tid],

            # adaptive priors - non-centered
            transpars> matrix[actor,4]:alpha <-
                    compose_noncentered( sigma_actor , L_Rho_actor , z_actor ),
            transpars> matrix[block_id,4]:beta <-
                    compose_noncentered( sigma_block , L_Rho_block , z_block ),
            matrix[4,actor]:z_actor ~ normal( 0 , 1 ),
            matrix[4,block_id]:z_block ~ normal( 0 , 1 ),

            # fixed priors
            g[tid] ~ normal(0,1),
            vector[4]:sigma_actor ~ dexp(1),
            cholesky_factor_corr[4]:L_Rho_actor ~ lkj_corr_cholesky( 2 ),
            vector[4]:sigma_block ~ dexp(1),
            cholesky_factor_corr[4]:L_Rho_block ~ lkj_corr_cholesky( 2 ),

            # compute ordinary correlation matrixes from Cholesky factors
            gq> matrix[4,4]:Rho_actor <<- Chol_to_Corr(L_Rho_actor),
            gq> matrix[4,4]:Rho_block <<- Chol_to_Corr(L_Rho_block)
        ) , data=dat , chains=4 , cores=4 , log_lik=TRUE )
    )
    expect_equivalent( dim(precis(m14.3,3)) , c(180,6) )
})

## R code 14.23
set.seed(73)
N <- 500
U_sim <- rnorm( N )
Q_sim <- sample( 1:4 , size=N , replace=TRUE )
E_sim <- rnorm( N , U_sim + Q_sim )
W_sim <- rnorm( N , U_sim + 0*E_sim )
dat_sim <- list(
    W=standardize(W_sim) ,
    E=standardize(E_sim) ,
    Q=standardize(Q_sim) )

## R code 14.24
test_that("R code 14.24",{
    m14.4 <- ulam(
        alist(
            W ~ dnorm( mu , sigma ),
            mu <- aW + bEW*E,
            aW ~ dnorm( 0 , 0.2 ),
            bEW ~ dnorm( 0 , 0.5 ),
            sigma ~ dexp( 1 )
        ) , data=dat_sim , chains=4 , cores=4 )
    expect_equivalent( dim(precis(m14.4,3)) , c(3,6) )
})

## R code 14.25
test_that("R code 14.25",{
    m14.5 <- ulam(
        alist(
            W ~ dnorm( mu , sigma ),
            mu <- aW + bEW*E + bQW*Q,
            aW ~ dnorm( 0 , 0.2 ),
            bEW ~ dnorm( 0 , 0.5 ),
            bQW ~ dnorm( 0 , 0.5 ),
            sigma ~ dexp( 1 )
        ) , data=dat_sim , chains=4 , cores=4 )
    expect_equivalent( dim(precis(m14.5,3)) , c(4,6) )
})

## R code 14.26
test_that("R code 14.26",{
    expect_warning(
    m14.6 <- ulam(
        alist(
            c(W,E) ~ multi_normal( c(muW,muE) , Rho , Sigma ),
            muW <- aW + bEW*E,
            muE <- aE + bQE*Q,
            c(aW,aE) ~ normal( 0 , 0.2 ),
            c(bEW,bQE) ~ normal( 0 , 0.5 ),
            Rho ~ lkj_corr( 2 ),
            Sigma ~ exponential( 1 )
        ), data=dat_sim , chains=4 , cores=4 )
    )
    expect_equivalent( dim(precis(m14.6,3)) , c(10,6) )
})

## R code 14.30
library(rethinking)
data(KosterLeckie)

## R code 14.31
kl_data <- list(
    N = nrow(kl_dyads),
    N_households = max(kl_dyads$hidB),
    did = kl_dyads$did,
    hidA = kl_dyads$hidA,
    hidB = kl_dyads$hidB,
    giftsAB = kl_dyads$giftsAB,
    giftsBA = kl_dyads$giftsBA
)

test_that("R code 14.31",{
    expect_warning(
    m14.7 <- ulam(
        alist(
            giftsAB ~ poisson( lambdaAB ),
            giftsBA ~ poisson( lambdaBA ),
            log(lambdaAB) <- a + gr[hidA,1] + gr[hidB,2] + d[did,1] ,
            log(lambdaBA) <- a + gr[hidB,1] + gr[hidA,2] + d[did,2] ,
            a ~ normal(0,1),

           ## gr matrix of varying effects
            vector[2]:gr[N_households] ~ multi_normal(0,Rho_gr,sigma_gr),
            Rho_gr ~ lkj_corr(4),
            sigma_gr ~ exponential(1),

           ## dyad effects
            transpars> matrix[N,2]:d <-
                    compose_noncentered( rep_vector(sigma_d,2) , L_Rho_d , z ),
            matrix[2,N]:z ~ normal( 0 , 1 ),
            cholesky_factor_corr[2]:L_Rho_d ~ lkj_corr_cholesky( 8 ),
            sigma_d ~ exponential(1),

           ## compute correlation matrix for dyads
            gq> matrix[2,2]:Rho_d <<- Chol_to_Corr( L_Rho_d )
        ), data=kl_data , chains=4 , cores=4 , iter=2000 )
    )
    expect_equivalent( dim(precis(m14.7,3)) , c(1266,6) )
})

## R code 14.37
# load the distance matrix
library(rethinking)
data(islandsDistMatrix)

# display (measured in thousands of km)
Dmat <- islandsDistMatrix
colnames(Dmat) <- c("Ml","Ti","SC","Ya","Fi","Tr","Ch","Mn","To","Ha")

## R code 14.39
data(Kline2) # load the ordinary data, now with coordinates
d <- Kline2
d$society <- 1:10 # index observations

dat_list <- list(
    T = d$total_tools,
    P = d$population,
    society = d$society,
    Dmat=islandsDistMatrix )

test_that("R code 14.39",{
    m14.8 <- ulam(
        alist(
            T ~ dpois(lambda),
            lambda <- (a*P^b/g)*exp(k[society]),
            vector[10]:k ~ multi_normal( 0 , SIGMA ),
            matrix[10,10]:SIGMA <- cov_GPL2( Dmat , etasq , rhosq , 0.01 ),
            c(a,b,g) ~ dexp( 1 ),
            etasq ~ dexp( 2 ),
            rhosq ~ dexp( 0.5 )
        ), data=dat_list , chains=4 , cores=4 , iter=2000 )
    expect_equivalent( dim(precis(m14.8,3)) , c(15,6) )
})

## R code 14.46
test_that("R code 14.46",{
    expect_warning(
    m14.8nc <- ulam(
        alist(
            T ~ dpois(lambda),
            lambda <- (a*P^b/g)*exp(k[society]),

            # non-centered Gaussian Process prior
            transpars> vector[10]: k <<- L_SIGMA * z,
            vector[10]: z ~ normal( 0 , 1 ),
            transpars> matrix[10,10]: L_SIGMA <<- cholesky_decompose( SIGMA ),
            transpars> matrix[10,10]: SIGMA <- cov_GPL2( Dmat , etasq , rhosq , 0.01 ),

            c(a,b,g) ~ dexp( 1 ),
            etasq ~ dexp( 2 ),
            rhosq ~ dexp( 0.5 )
        ), data=dat_list , chains=4 , cores=4 , iter=2000 )
    )
    expect_equivalent( dim(precis(m14.8nc,3)) , c(225,6) )
})

## R code 14.47
library(rethinking)
data(Primates301)
data(Primates301_nex)

# plot it using ape package - install.packages('ape') if needed
library(ape)

## R code 14.48
d <- Primates301
d$name <- as.character(d$name)
dstan <- d[ complete.cases( d$group_size , d$body , d$brain ) , ]
spp_obs <- dstan$name

## R code 14.49
dat_list <- list(
    N_spp = nrow(dstan),
    M = standardize(log(dstan$body)),
    B = standardize(log(dstan$brain)),
    G = standardize(log(dstan$group_size)),
    Imat = diag(nrow(dstan)) )

test_that("R code 14.49",{
    m14.9 <- ulam(
        alist(
            B ~ multi_normal( mu , SIGMA ),
            mu <- a + bM*M + bG*G,
            matrix[N_spp,N_spp]: SIGMA <- Imat * sigma_sq,
            a ~ normal( 0 , 1 ),
            c(bM,bG) ~ normal( 0 , 0.5 ),
            sigma_sq ~ exponential( 1 )
        ), data=dat_list , chains=4 , cores=4 )
    expect_equivalent( dim(precis(m14.9,3)) , c(4,6) )
})

## R code 14.50
library(ape)
tree_trimmed <- keep.tip( Primates301_nex, spp_obs )
Rbm <- corBrownian( phy=tree_trimmed )
V <- vcv(Rbm)
Dmat <- cophenetic( tree_trimmed )

## R code 14.51
# put species in right order
dat_list$V <- V[ spp_obs , spp_obs ]
# convert to correlation matrix
dat_list$R <- dat_list$V / max(V)

# Brownian motion model
test_that("R code 14.51",{
    m14.10 <- ulam(
        alist(
            B ~ multi_normal( mu , SIGMA ),
            mu <- a + bM*M + bG*G,
            matrix[N_spp,N_spp]: SIGMA <- R * sigma_sq,
            a ~ normal( 0 , 1 ),
            c(bM,bG) ~ normal( 0 , 0.5 ),
            sigma_sq ~ exponential( 1 )
        ), data=dat_list , chains=4 , cores=4 )
    expect_equivalent( dim(precis(m14.10,3)) , c(4,6) )
})

## R code 14.52
# add scaled and reordered distance matrix
dat_list$Dmat <- Dmat[ spp_obs , spp_obs ] / max(Dmat)

test_that("R code 14.52",{
    m14.11 <- ulam(
        alist(
            B ~ multi_normal( mu , SIGMA ),
            mu <- a + bM*M + bG*G,
            matrix[N_spp,N_spp]: SIGMA <- cov_GPL1( Dmat , etasq , rhosq , 0.01 ),
            a ~ normal(0,1),
            c(bM,bG) ~ normal(0,0.5),
            etasq ~ half_normal(1,0.25),
            rhosq ~ half_normal(3,0.25)
        ), data=dat_list , chains=4 , cores=4 )
    expect_equivalent( dim(precis(m14.11,3)) , c(5,6) )
})


