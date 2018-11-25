# simulation for Chapter 6 collider example

sim_happiness <- function( seed=1977 , N_years=1000 , max_age=65 , N_births=20 , aom=18 ) {
    set.seed(seed)
    H <- M <- A <- c()
    for ( t in 1:N_years ) {
        A <- A + 1 # age existing individuals
        A <- c( A , rep(1,N_births) ) # newborns
        H <- c( H , seq(from=-2,to=2,length.out=N_births) ) # sim happiness trait - never changes
        M <- c( M , rep(0,N_births) ) # not yet married
        # for each person over 17, chance get married
        for ( i in 1:length(A) ) {
            if ( A[i] >= aom & M[i]==0 ) {
                M[i] <- rbern(1,inv_logit(H[i]-4))
            }
        }
        # mortality
        deaths <- which( A > max_age )
        if ( length(deaths)>0 ) {
            A <- A[ -deaths ]
            H <- H[ -deaths ]
            M <- M[ -deaths ]
       }
    }
    d <- data.frame(age=A,married=M,happiness=H)
    return(d)
}
