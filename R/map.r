# MAP - maximum a posteriori
# uses Bolker's mle2 (bbmle) as the engine, just smuggling in priors as penalized likelihood

# format for priors:
# list( formula1 , formula2 , ... )
# formulas are in format of the likelihood formula, i.e.
# parameter ~ dfunction( pars )
# e.g.:
# a ~ dnorm(0,1)
# s ~ dunif(0,100)

# map() then builds the likelihood formula and prior formulas into a unified enclosure that returns penalized log-likelihood. optim() does the rest.

map <- function( formula , data , start , prior , ... ) {

    require(bbmle)
    
    # check parameters
    
    if ( missing( formula ) | missing( data ) | missing( start ) ) {
        stop( "Must specify formula, data, and start." )
    }
    if ( missing( prior ) ) {
        # okay, but mark as NULL
        prior <- NULL
    }
    
    # check that none of the declared priors are for a parameter not in start list
    
    #################
    # build enclosure
    
}

# testing examples
# data(cars)
# fit <- map( dist ~ dnorm( mean=a+b*speed , sd=sigma ) , data=cars , start=list( a=mean(cars$dist) , b=0 , sigma=sd(cars$dist) ) , prior=list( b ~ dnorm(0,1) , sigma ~ dunif(0,100) ) )