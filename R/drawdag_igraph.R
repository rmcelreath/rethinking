# my DAG drawing function that uses igraph library

graphdag <- function( x , 
    layout=layout_nicely ,  
    vertex.size=30 , 
    vertex.label.family="sans" ,
    edge.curved=0 , 
    edge.color="black" , 
    edge.width=1.5,
    edge.arrow.width=1.2 , 
    edge.arrow.size=0.7 ,
    interact=FALSE ,
    ... ) {
    require(igraph)
    x <- as.dagitty(x)
    # extract edges and nodes
    the_edges <- dagitty::edges(x)
    the_vars <- unique( c( levels(the_edges[[1]]) , levels(the_edges[[2]]) ) )
    # translate edges list to adjacency matrix
    n <- length(the_vars)
    adjmat <- adjmat_from_dag( x )
    # existing coordinates?
    coords <- coordinates( x )
    if ( all( !is.na(coords$x) ) ) {
        # coordinates complete - convert to matrix
        # y axis in dagitty is flipped
        layout <- matrix( NA , nrow=n , ncol=2 )
        for ( i in 1:n ) {
            xx <- coords$x[[ the_vars[i] ]]
            yy <- (-1)*coords$y[[ the_vars[i] ]]
            layout[i,] <- c(xx,yy)
        }
    }
    # draw
    mgraph <- graph_from_adjacency_matrix( adjmat , mode="directed" )
    the_shapes <- rep("none",n)
    # any latent vars?
    unobs_vars <- dagitty::latents(x)
    if ( length(unobs_vars)>0 ) {
        for ( uv in unobs_vars ) the_shapes[ which(the_vars==uv) ] <- "circle"
    }
    if ( interact==FALSE ) {
        if ( class(layout)!="matrix" )
            the_layout <- do.call( layout , list(mgraph) )
        else
            the_layout <- layout
        plot( mgraph , 
            vertex.color = "white" , 
            vertex.size = vertex.size , 
            vertex.shape = the_shapes ,
            vertex.label.family = vertex.label.family ,
            edge.arrow.size = edge.arrow.size , 
            edge.arrow.width = edge.arrow.width ,
            edge.curved = edge.curved , 
            edge.color = edge.color ,
            edge.width = edge.width ,
            vertex.label = the_vars , 
            vertex.label.color = "black" ,
            seed = 1 , layout=the_layout , ... )
    } else {
        tkplot( mgraph )
    }
    rownames(the_layout) <- the_vars
    colnames(the_layout) <- c("x-coord","y-coord")
    return(invisible(the_layout))
}

# take dagitty dag and translate to adjacency matrix
adjmat_from_dag <- function( x ) {
    the_edges <- dagitty::edges(x)
    the_vars <- unique( c( levels(the_edges[[1]]) , levels(the_edges[[2]]) ) )
    # translate edges list to adjacency matrix
    n <- length(the_vars)
    adjmat <- matrix( 0 , nrow=n , ncol=n )
    rownames(adjmat) <- the_vars
    colnames(adjmat) <- the_vars
    for ( i in 1:nrow(the_edges) ) {
        node1 <- as.character( the_edges[i,1] )
        node2 <- as.character( the_edges[i,2] )
        tie_type <- as.character( the_edges[i,3] )
        if ( tie_type=="->" ) {
            # make directed tie
            adjmat[ node1 , node2 ] <- 1
        }
    }
    return(adjmat)
}

# function to interactively draw layout for given dag
sketchdag <- function( x , cleanup=1 , plot=TRUE , rescale=FALSE , grid=0.2 , ... ) {
    the_edges <- dagitty::edges(x)
    the_vars <- unique( c( levels(the_edges[[1]]) , levels(the_edges[[2]]) ) )
    n <- length(the_vars)
    adjmat <- adjmat_from_dag( x )
    print(the_vars)
    plot( NULL , xlim=c(-1,1) , ylim=c(-1,1) , type='n' , axes=FALSE , ann=FALSE )
    if ( grid != FALSE ) {
        # draw alignment grid
        cseq <- seq( from=-1 , to=1 , by=grid )
        for ( xc in cseq ) abline( h=xc , lwd=0.5 )
        for ( yc in cseq ) abline( v=yc , lwd=0.5 )
    }
    pts <- matrix( NA , nrow=n , ncol=2 )
    for ( i in 1:n ) {
        pt <- locator( 1 , type="p" , pch=the_vars[i] )
        pts[i,1] <- pt$x
        pts[i,2] <- pt$y
    }
    rownames(pts) <- the_vars
    colnames(pts) <- c("x-coord","y-coord")
    if ( cleanup != FALSE ) pts <- round( pts , cleanup )
    if ( plot==TRUE )
        return(graphdag(x,layout=pts,...))
    else
        return(pts)
}

# test code
if (FALSE) {

    #Other graph layouts: add_layout_, component_wise, layout_as_bipartite, layout_as_star, layout_as_tree, layout_in_circle, layout_nicely, layout_on_grid, layout_on_sphere, layout_randomly, layout_with_dh, layout_with_fr, layout_with_gem, layout_with_graphopt, layout_with_lgl, layout_with_mds, layout_with_sugiyama, layout_, merge_coords, norm_coords, normalize

    library(rethinking)

    exdag <- dagitty( "dag {
        U [unobserved]
        X [exposure]
        Y [outcome]
        Z -> X -> Y
        X <- U -> Y
    }")

    coordinates(exdag) <- list(x=c(U=1,X=1,Y=2,Z=0),y=c(U=0,X=1,Y=1,Z=1))
    graphdag(exdag,margin=-0.25,edge.curve=0.1 )

    l <- graphdag( exdag , layout=layout_in_circle )
    l <- graphdag( exdag , layout=layout_nicely )
    graphdag( exdag , layout=l )

    exdag2 <- dagitty( 'dag {
        G [exposure,unobserved]
        P [outcome]
        "G*" <- G -> P -> W -> RG
        RG -> "G*"
    }')

    graphdag( exdag2 , edge.curve=0.2 )

    l <- sketchdag( exdag2 )
    drawdag( exdag2 , layout=round(l,2) )

    sketchdag( exdag2 , asp=0.8 )

}
