# IC class and show method

setClass( "ICnumeric" , contains="numeric" )

setMethod( "show" , "ICnumeric" , function(object) {
    cat( round( object , 1 ) )
    if ( !is.null( attr(object,"se") ) ) {
        cat( concat( " ±" , round( attr(object,"se") , 1 ) ) )
    }
    cat("\n")
} )
