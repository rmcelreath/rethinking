# lm link and sim methods
# later to add glm

setMethod("link", "lm",
function( fit , data , n=1000 , post , refresh=0.1 , replace=list() , flatten=TRUE , ... ) {

    # extract samples --- will need for inline parameters e.g. sigma in likelihood
    if ( missing(post) ) 
        post <- extract.samples(fit,n=n)
    else {
        n <- dim(post[[1]])[1]
        if ( is.null(n) ) n <- length(post[[1]])
    }
    
    # replace with any elements of replace list
    if ( length( replace ) > 0 ) {
        for ( i in 1:length(replace) ) {
            post[[ names(replace)[i] ]] <- replace[[i]]
        }
    }
    
    # check for new data frame
    # if there, replace fit$model with it
    # BUT must check for outcome variable and insert
    # otherwise predict() won't work later
    if ( !missing(data) ) {
        model_orig <- fit$model
        out_var <- names(model_orig)[1] # first name should be outcome
        data <- cbind( 0 , data )
        names(data)[1] <- out_var
        fit$model <- data
    }
    
    # compute mu at each sample from posterior
    # can hijack predict.lm() to do this by inserting samples into object$coefficients
    m_temp <- fit
    out_var <- fit$model[[1]] # should be outcome variable
    n_coefs <- length(coef(fit))
    mu <- sapply( 1:n , 
        function(i) {
            for ( j in 1:n_coefs ) m_temp$coefficients[[j]] <- post[i,j]             
            as.numeric(predict(m_temp)) # as.numeric strips names
        } )
    
    return(t(mu)) # need to t(), so rows are samples and cols are obs
} )

setMethod("sim", "lm",
function( fit , data , n=1000 , post , ll=FALSE , refresh=0.1 , ... ) {

    # extract samples --- will need for inline parameters e.g. sigma in likelihood
    if ( missing(post) ) 
        post <- extract.samples(fit,n=n)
    else {
        n <- dim(post[[1]])[1]
        if ( is.null(n) ) n <- length(post[[1]])
    }
    
    # check for new data frame
    # if there, replace fit$model with it
    # BUT must check for outcome variable and insert
    # otherwise predict() won't work later
    if ( !missing(data) ) {
        model_orig <- fit$model
        out_var <- names(model_orig)[1] # first name should be outcome
        data <- cbind( 0 , data )
        names(data)[1] <- out_var
        fit$model <- data
    }
    
    # get sigma
    sigma <- sd(resid(fit))
    
    # draw samples obs at each sample from posterior
    # can hijack predict.lm() to do this by inserting samples into object$coefficients
    m_temp <- fit
    out_var <- fit$model[[1]] # should be outcome variable
    n_coefs <- length(coef(fit))
    y <- sapply( 1:n , 
        function(i) {
            for ( j in 1:n_coefs ) m_temp$coefficients[[j]] <- post[i,j]             
            mu <- as.numeric(predict(m_temp)) # as.numeric strips names
            return( rnorm(length(mu),mu,sigma) )
        } )
    
    return(t(y)) # need to t(), so rows are samples and cols are obs
} )
