\name{Primates301}
\alias{Primates301}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Primate life history and social learning data}
\description{
  Life history data, social learning data, and phylogenetic distance matrix for 301 primate species. These data were assembled by and analyzed in Street et al 2017 (see references). 
}
\usage{
data(Primates301)
data(Primates301_distance_matrix)
data(Primates301_vcov_matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
}
\format{
    \code{Primates301} is a data.table with elements: 
    \enumerate{
        \item name: Full taxonomic name of species
        \item genus : Genus of species
        \item species : Species name within genus
        \item subspecies : Sub-species designation, if any
        \item spp_id : Unique ID for species
        \item genus_id : Unique ID for genus
        \item social_learning : Count of mentions of social learning in literature
        \item research_effort : Size of literature on species
        \item brain : Brain volume (endocranial volume) in cubic centimeters
        \item body : Body mass in grams
        \item group_size : Average social group size
        \item gestation : Length of gestation (days)
        \item weaning : At at weaning (days)
        \item longevity : Maximum lifespan (months)
        \item sex_maturity : Age of sexual maturity (days)
        \item maternal_investment : Period of maternal investment (days) = gestation + weaning
    }
    \code{Primates301_distance_matrix} is a matrix with species on the margins and phylogenetic distances in the cells.
    \code{Primates301_vcov_matrix} is a matrix with species on the margins and variances-covariances in the cells.
}
\value{
}
\references{
Street SE, Navarrete AF, Reader SM, Laland KN (2017) Coevolution of cultural intelligence, extended life history, sociality, and brain size in primates. PNAS https://doi.org/10.1073/pnas.1620734114

Arnold C, Matthews LJ, Nunn CL (2010) The 10kTrees Website: A New Online Resource for Primate Phylogeny. Evol Anthropol 19(3):114-118.

Reader SM, Hager Y, Laland KN (2011) The evolution of primate general and cultural intelligence. Philos Trans R Soc B-Biological Sci 366(1567):1017-1027.

Isler K, et al. (2008) Endocranial volumes of primate species: scaling analyses using a comprehensive and reliable data set. J Hum Evol 55(6):967-978.

Jones, Kate E, et al. (2009) PanTHERIA: a species-level database of life history, ecology, and geography of extant and recently extinct mammals. Ecology 90:2649.
}
\seealso{}
\examples{
data(Primates301)
plot( log(brain) ~ log(body) , data=Primates301 )

data(Primates301_distance_matrix)
image(Primates301_distance_matrix)

# Gaussian process phylogenetic regression
# prep variables
d <- Primates301
d$name <- as.character(d$name)
dstan <- d[ complete.cases( d$social_learning, d$research_effort , d$body , d$brain ) , ]
# prune distance matrix to spp in dstan
spp_obs <- dstan$name
y <- Primates301_distance_matrix
y2 <- y[ spp_obs , spp_obs ]
# cbind( sort(spp_obs) , sort(colnames(y2)) )
# scale distances
y3 <- y2/max(y2)

mP301GP <- ulam(
    alist(
        social_learning ~ poisson( lambda ),
        log(lambda) <- a + g[spp_id] + b_ef*log_research_effort + b_body*log_body + b_eq*log_brain,
        a ~ normal(0,1),
        vector[N_spp]: g ~ multi_normal( 0 , SIGMA ),
        matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat , etasq , rhosq , 0.01 ),
        b_body ~ normal(0,1),
        b_eq ~ normal(0,1),
        b_ef ~ normal(1,1),
        etasq ~ exponential(1),
        rhosq ~ exponential(1)
    ),
    data=list(
        N_spp = nrow(dstan),
        social_learning = dstan$social_learning,
        spp_id = 1:nrow(dstan),
        log_research_effort = log(dstan$research_effort),
        log_body = log(dstan$body),
        log_brain = log(dstan$brain),
        Dmat = y3
    ) , 
    control=list(max_treedepth=15,adapt_delta=0.95) , 
    sample=FALSE , iter=400 )

# non-centered, Cholesky form
mP301GPnc <- ulam(
    alist(
        social_learning ~ poisson( lambda ),
        log(lambda) <- a + g[spp_id] + b_ef*log_research_effort + b_body*log_body + b_eq*log_brain,
        a ~ normal(0,1),
        vector[N_spp]: g <<- L_SIGMA * eta,
        vector[N_spp]: eta ~ normal( 0 , 1 ),
        matrix[N_spp,N_spp]: L_SIGMA <<- cholesky_decompose( SIGMA ),
        matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat , etasq , rhosq , 0.01 ),
        b_body ~ normal(0,1),
        b_eq ~ normal(0,1),
        b_ef ~ normal(1,1),
        etasq ~ exponential(1),
        rhosq ~ exponential(1)
    ),
    data=list(
        N_spp = nrow(dstan),
        social_learning = dstan$social_learning,
        spp_id = 1:nrow(dstan),
        log_research_effort = log(dstan$research_effort),
        log_body = log(dstan$body),
        log_brain = log(dstan$brain),
        Dmat = y3
    ) , 
    control=list(max_treedepth=15,adapt_delta=0.95) , 
    sample=FALSE , iter=400 )

# Pagel's lambda approach --- Not endorsed!
# This is of historical interest only
data(Primates301_vcov_matrix)
vcov_thin <- Primates301_vcov_matrix[ spp_obs , spp_obs ]
mP301L <- ulam(
    alist(
        social_learning ~ poisson( lambda ),
        log(lambda) <- a + g[spp_id] + b_ef*log_research_effort + b_body*log_body + b_eq*log_brain,
        a ~ normal(0,1),
        vector[N_spp]: g <<- L_SIGMA * eta,
        vector[N_spp]: eta ~ normal( 0 , 1 ),
        matrix[N_spp,N_spp]: L_SIGMA <<- cholesky_decompose( SIGMA ),
        matrix[N_spp,N_spp]: SIGMA <- cov_Pagel( SIGMA_raw , Plambda ),
        b_body ~ normal(0,1),
        b_eq ~ normal(0,1),
        b_ef ~ normal(1,1),
        Plambda ~ beta(2,2)
    ),
    data=list(
        N_spp = nrow(dstan),
        social_learning = dstan$social_learning,
        spp_id = 1:nrow(dstan),
        log_research_effort = log(dstan$research_effort),
        log_body = log(dstan$body),
        log_brain = log(dstan$brain),
        SIGMA_raw = vcov_thin
    ) , 
    control=list(max_treedepth=15,adapt_delta=0.95) , 
    sample=TRUE , iter=400 )

# centered version --- seems to mix better
mP301L2 <- ulam(
    alist(
        social_learning ~ poisson( lambda ),
        log(lambda) <- a + g[spp_id] + b_ef*log_research_effort + b_body*log_body + b_eq*log_brain,
        a ~ normal(0,1),
        vector[N_spp]: g ~ multi_normal( 0 , SIGMA ),
        matrix[N_spp,N_spp]: SIGMA <- cov_Pagel( SIGMA_raw , Plambda ),
        b_body ~ normal(0,1),
        b_eq ~ normal(0,1),
        b_ef ~ normal(1,1),
        Plambda ~ beta(2,2)
    ),
    data=list(
        N_spp = nrow(dstan),
        social_learning = dstan$social_learning,
        spp_id = 1:nrow(dstan),
        log_research_effort = log(dstan$research_effort),
        log_body = log(dstan$body),
        log_brain = log(dstan$brain),
        SIGMA_raw = vcov_thin
    ) , 
    control=list(max_treedepth=15,adapt_delta=0.95) , 
    sample=TRUE , iter=400 )

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

