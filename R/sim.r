# sim -- simulates data using formulas similar to map and map2stan

sim <- function( flist , values , debug=FALSE , ... ) {
    
    ########################################
    # check arguments
    if ( missing(flist) ) stop( "Formula required." )
    if ( class(flist) != "list" ) {
        if ( class(flist)=="formula" ) {
            flist <- list(flist)
        } else {
            stop( "Formula or list of formulas required." )
        }
    }
    if ( missing(values) ) stop( "'values' list required." )
    
}

# EXAMPLES
if ( FALSE ) {

f <- list(
    y ~ dnorm( mu , sigma )
)

} #EXAMPLES