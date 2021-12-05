# animation functions

ani.saveqz <- function (list , dpi=150 , location="frames/frame_" , prefix=1000 ) {
    require(animation)   
    if (missing(list)) 
        list = animation:::.ani.env$.images
    lapply(1:length(list), function(x) {
        dev.hold()
        replayPlot( list[[x]] )
        ani.pause()
        quartz.save( concat(location,prefix+x,".png") , type = "png", device = dev.cur(), dpi = dpi )
    })
    invisible(NULL)
}


anim_prior_predictive <- function(
    prior, # named list samples from prior
    linkf, # function (of x) to display 
    n=20, # number of samples to show
    xlim=c(-2,2), # range of x to show
    ylim=c(0,1), # y range
    n_pts=100, # number of points in each curve
    n_steps=10, # number of transition frames between subsequent samples
    n_to_show=1, # number of samples to show at once
    cols=c(2,4,5,6,3),
    alpha=0.7,
    xlab="x value",
    ylab="probability",
    lwd=5,
    pf, # optional function to control setting up plot and axes
    vel0=1, accel=-1, # spring parameters
    start_from, # optional initial parameter values (for transition from previous seq)
    do_reset=TRUE # reset animation stack
) {

    # function for 1D path acceleration - distance traveled at time t
    faccel <- function(t,x0=0,v0=vel0,a=accel) x0 + v0*t + (1/2)*a*t^2
    # interpolate between parameter values with deceleration
    interp <- function(p1,p2) {
        path <- sapply( seq(from=0,to=1,length=n_steps/2+1) , faccel )
        path_seq <- normalize(path)
        nn <- length(path_seq)
        path_seq <- c( 0.5*(1-path_seq)[nn:2] , 0.5*path_seq[1:(nn-1)] + 0.5 )
        #the_seq <- seq( from=s1 , to=s2 , length=n_steps )
        path_seq <- path_seq*p2 + (1-path_seq)*p1
        return(path_seq)
    }

    if ( missing(pf) ) 
        pf <- function(samp,step) {
            plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,bty="n")
            #axis(2,at=c(0,0.5,1),labels=c(0,0.5,1))
        }

    # setup colors
    if ( length(cols) < n_to_show ) cols <- rep(cols,length.out=n_to_show)
    cols <- sapply( cols , function(f) col.alpha(f,alpha) )

    # x values
    x <- seq( from=xlim[1] , to=xlim[2] , length.out=n_pts )

    require(animation)
    ani.record(reset = do_reset)  # clear history before recording

    # set up data structure for interpolation
    prior_blend <- list()
    for ( i in 1:length(prior) ) 
        prior_blend[[ names(prior)[i] ]] <- matrix(NA,nrow=n_to_show,ncol=n_steps)

    # check for start_from
    if ( !missing(start_from) )
        if ( !is.null(start_from) ) {
            # insert start_from into prior
            for ( j in names(prior) )
                for ( k in 1:n_to_show ) {
                    prior[[j]][ 1 + (k-1)*n ] <- start_from[[j]][ k , n_steps ]
                }
        }

    # animation loop
    for ( i in 1:(n-1) ) {

        # interpolate
        for ( j in names(prior) )
            for ( k in 1:n_to_show ) {
                s1 <- prior[[j]][ i + (k-1)*n ]
                s2 <- prior[[j]][ i+1 + (k-1)*n ]
                # prior_blend[[j]][k,] <- seq( from=s1 , to=s2 , length=n_steps )
                prior_blend[[j]][k,] <- interp( s1 , s2 )
            }

        for ( j in 1:n_steps ) {

            pf(i,j)

            for ( k in 1:n_to_show ) {

                prior_now <- list()
                for ( kk in names(prior) ) prior_now[[kk]] <- prior_blend[[kk]][k,j]
                prior_now[['x']] <- x

                y <- eval( body(linkf) , env=prior_now )
                lines( x , y , lwd=lwd , col=cols[k] )
            }
            
            ani.record()  # record the current frame

        }
        
    }#i
    return(prior_blend)
}


#### DAG STUFF ####


# function to draw compound effect dot at particular coordinate
dagfx_dot <- function(x,y,k,r=1,rr=1,angle1=pi/2,lwd=10,angle_gap=pi/5,pch=16) {
    # k should be vector of colors
    require(plotrix)
    n <- length(k)
    d <- rep(NA,n)
    d[1] <- r
    #if ( n > 1 ) for ( i in 2:n ) d[i] <- d[i-1] + rr
    if ( n > 1 ) for ( i in 2:n ) d[i] <- d[1] + rr
    angle_inc <- 2*pi / (n-1)
    #angle1 <- pi/2
    #angle_gap <- pi/8
    for ( i in n:1 ) {
        if ( i==1 ) {
            # draw center circle
            points( x , y , pch=16 , cex=d[i]*1.2 , col="black" )
            points( x , y , pch=pch , cex=d[i] , col=k[i] )
        } else {
            # draw arc
            draw.arc( x , y , radius=d[i]/80 , col=k[i] , lwd=lwd , angle1=angle1 + angle_gap , angle2=angle1 + angle_inc - angle_gap )
            angle1 <- angle1 + angle_inc
        }
    }
}
# plot( NULL , xlim=c(-1,1) , ylim=c(-1,1) ); dagfx_dot(0,0,2:4 , r=2 , lwd=8 )

# master function
# starts animation and tracks compound effects through FRONT DOORS
dagfx_anim_forward <- function( the_dag , Y , X , n_frames , n_loops=3 , path_frames=20 , radius=2 , color_list=c(2,4,3,5,6,7) , fade_rate=0.9 , col="white" , bg="black" , xlim , ylim , buffer=c(0,0.4,0,0) , offset_stack=0.12 , arc_lwd=10 , angle_gap=pi/5 , angle1_inc=pi/path_frames , fix="" , skadoosh=TRUE , accel=-0.7 , vel0=1 , force_on="" , scm , ... ) {

    # function to extract colors for each dot
    dot_palette <- function( parent ) {
        # get parent color
        the_col <- color_list[ which( uvars == parent ) ]
        # append colors of vars that touched parent
        if ( any(touched[ parent , ]>0) ) {
            j <- which( touched[ parent , ] > 0 )
            for ( k in j ) {
                a <- touched[ parent , k ]
                the_col <- c( the_col , col.alpha( color_list[ k ] , a ) )
            }
        }
        return(the_col)
    }

    # set frame count
    if ( missing(n_frames) & !missing(n_loops) ) {
        n_frames <- n_loops * path_frames
    }

    # get direct animation paths
    the_edges <- edges(the_dag)
    the_edges <- the_edges[ the_edges[,3]=="->" , 1:3 ]
    the_edges[,3] <- the_edges$v

    # calculate plot dimensions so there is space for accumulating on outcome
    coords <- coordinates(the_dag)
    labels <- names(coords$x)
    wx <- sapply( paste0("mm",labels), 
        function(s) strwidth(s,units="inches") )
    wy <- sapply( paste0("\n",labels), 
        function(s) strheight(s,units="inches") )
    if ( missing(xlim) )
        xlim <- c( min(coords$x-wx/2) , max(coords$x+wx/2) ) + c( -buffer[1] , buffer[2] )
    if ( missing(ylim) )
        ylim <- c( -max(coords$y+wy/2) , -min(coords$y-wy/2) ) + c( -buffer[3] , buffer[4] )

    # draw and collect path coordinates
    par(bg = bg)
    coo_arrows <- drawdag( the_dag , col_arrow=col , col_segment=col , col_labels=col , xlim=xlim, ylim=ylim )
    coo_nodes <- coordinates( the_dag )

    # check fix (conditioning) list for colliders
    # if any colliders in fix list, need to add new paths for them
    # idea is to "reflect" parents of collider to other side of the collider
    # so e.g. A -> Z <- B
    # when Z is fixed, then A gets a dose of backdoor B and the reverse
    is_collider <- function(z) {
        # if z is endpoint for more than one edge, then collider
        if ( length( which( the_edges$w==z ) ) > 1 ) {
            #print( concat( z , " is a collider" ) )
            return(TRUE)
        }
        return(FALSE)
    }
    if ( !missing(fix) ) {
        for ( v in fix )
            if ( is_collider(v) ) {
                # turn off normal outbound (frontdoor) paths for v
                # we can do this instead by not making them "touch"
                ew <- which( the_edges$v==v )
                #the_edges[ ew , 3 ] <- NA
                # add edges to reflect v back on appropriate parents
                ew <- which( the_edges$w==v )
                collider_parents <- the_edges$v[ew]
                for ( w in collider_parents ) {
                    for ( ww in collider_parents )
                        if ( ww != w ) {
                            # add ghost path from collider v to parent w, which transmits opposite parent ww
                            cur_edg <- nrow(the_edges) + 1
                            the_edges[ cur_edg , ] <- c( v , w , ww )
                            # also copy coordinates into coo_arrows
                            # can copy path w -> v and then reverse it
                            old_arrow_i <- which( the_edges$v==w & the_edges$w==v )
                            new_arrow <- coo_arrows[old_arrow_i, c( 3 , 4 , 1 , 2 , 5 ) ]
                            coo_arrows[ cur_edg , ] <- new_arrow
                        }
                }#w
            }
        print( the_edges )
        print( coo_arrows )
    }

    # need list of which variables have "touched" others so we can draw the dots right
    uvars <- unique( c( unique(the_edges[,1]) , unique(the_edges[,2]) ) )
    touched <- matrix(0,nrow=length(uvars),ncol=length(uvars))
    colnames(touched) <- uvars
    rownames(touched) <- uvars

    # data frame, for when scm present
    dat <- matrix(NA,nrow=n_loops+2,ncol=length(uvars))
    colnames(dat) <- uvars
    dat_pos <- 0
    dat_sampler <- "" # record which parent updates outcome first

    # sample variables
    for ( i in 1:nrow(dat) ) {
        sv <- rep(NA,length(uvars))
        names(sv) <- uvars
        for ( u in scm ) {
            sv[u$var] <- u$f(sv)
            dat[i,u$var] <- sv[u$var]
        }
    }#i
    print(dat)

    # organize colors
    if ( class(color_list)=="list" ) {
        # named list, so sort in order of uvars
        color_list <- unlist(color_list)[uvars]
    }

    # positions of each effect
    pos_inc <- 1/path_frames
    #pos <- rep( -pos_inc , nrow(the_edges) )
    pos <- rep( 0 , nrow(the_edges) ) # init position

    # matrix with same dims as touched, used to track things that arrive at Y
    Ystack <- touched

    # path sequence
    # function for 1D path acceleration - distance traveled at time t
    fa <- function(t,x0=0,v0=vel0,a=accel) x0 + v0*t + (1/2)*a*t^2
    path_seq <- seq( from=0 , to=1 , length=path_frames )
    if ( skadoosh==TRUE ) {
        x <- sapply( path_seq , fa )
        path_seq <- normalize(x)
    }
    path_seq <- rep( path_seq , length=n_frames )

    # animation loop
    angle1 <- pi/2
    #angle1_inc <- pi/path_frames
    ani.record(reset = TRUE)
    for ( f in 1:n_frames ) {

        # draw the dag again
        par(bg = bg)
        drawdag( the_dag , col_arrow=col , col_segment=col , col_labels=col , xlim=xlim, ylim=ylim )

        # draw boxes on fixed variables
        if ( !missing(fix) ) {
            for ( v in fix ) {
                points( coo_nodes$x[v] , -coo_nodes$y[v] , pch=0 , cex=3 , col=col )
            }
        }

        # fade touched
        touched <- touched * fade_rate

        # for each effect, update dot position and draw
        for ( i in 1:nrow(the_edges) ) {

            parent <- the_edges$v[i]
            child <- the_edges$w[i]

            #if ( !missing(fix) ) {
            #    if ( parent %in% fix & !is_collider(parent) ) {
            #        # skip this (conditioned mediator or confound)
            #        next
            #    }
            #} 

            turned_off <- parent %in% fix # won't touch child
            if ( ( 
                length(parents(the_dag,parent))==0 | 
                any(touched[parent,]>0) | 
                parent %in% force_on 
            ) ) {

                #pos[i] <- pos[i] + pos_inc
                pos[i] <- pos[i] + 1
                if ( path_seq[pos[i]] == 1 ) {
                    # path end!
                    # run into descendant and pass on our color
                    # w is descendant, v is parent
                    # need to check e column of the_edges to know which var does the touching
                    if ( !missing(fix) )
                        if ( !(the_edges$e[i] %in% fix) )
                            touched[ the_edges$w[i] , the_edges$e[i] ] <- 1
                    # reset position to origin
                    #pos[i] <- (-pos_inc)
                    # if path end is the Y (monitored outcome),
                    # then stack
                    if ( !missing(Y) ) {
                        if ( Y == child ) {
                            Ystack[parent,] <- touched[parent,]
                            Ystack[parent,parent] <- 1
                        
                            if ( !missing(scm) & ( dat_sampler==parent | dat_sampler=="" ) ) {
                                # update counter for number of rows to display in stack
                                dat_pos <- dat_pos + 1
                                dat_sampler <- parent
                            }#scm
                        }#Y child
                    }#Y end
                }#path end

                # interpolate position along path
                px <- (1-path_seq[pos[i]])*coo_arrows[i,1] + path_seq[pos[i]]*coo_arrows[i,3]
                py <- (1-path_seq[pos[i]])*coo_arrows[i,2] + path_seq[pos[i]]*coo_arrows[i,4]

                # get parent color
                the_col <- dot_palette( the_edges$e[i] )
                #the_col <- color_list[ which( uvars == parent ) ]
                # append colors of vars that touched parent
                #if ( any(touched[ parent , ]>0) ) {
                #    j <- which( touched[ parent , ] > 0 )
                #    for ( k in j ) {
                #        a <- touched[ parent , k ]
                #        the_col <- c( the_col , col.alpha( color_list[ k ] , a ) )
                #    }
                #}

                # draw dot
                dot_type <- 16
                dagfx_dot( px , -py , the_col , r=radius , angle1=angle1 , lwd=arc_lwd , angle_gap=angle_gap , pch=dot_type )

            }

        }#i

        # draw the Ystack
        if ( !missing(Y) ) {
            pos_stack <- 1 # position in visual stack
            for ( v in rownames(Ystack) ) {
                if ( v != Y ) {
                    if ( Ystack[v,v] == 1 ) {
                        # hit, so display
                        px <- coo_nodes$x[Y]
                        py <- coo_nodes$y[Y]
                        # get parent color
                        the_col <- dot_palette( v )
                        # draw dot
                        j <- ifelse( length(the_col) > 1 , 1 , 0 )
                        px <- px + offset_stack*pos_stack + 0.03*j
                        dagfx_dot( px , -py , the_col , r=radius , angle1=angle1 , lwd=arc_lwd  , angle_gap=angle_gap )
                        # draw label above
                        py <- -py + 0.12
                        text( px , py , v , col=col )
                        # update position
                        pos_stack <- pos_stack + 1 + 0.4*j
                    }
                }
            }#i
        }

        angle1 <- angle1 + angle1_inc

        # record frame
        ani.record()

    }#f

    # print(coo_arrows)
    # print(touched)
    # print(Ystack)

    return(invisible(dat))

}

#######
# DAG animations using backdoor paths

dagfx_anim_backward <- function( the_dag , Y , X , n_frames , n_loops=3 , path_frames=20 , radius=2 , color_list=c(2,4,3,5,6,7) , col="white" , bg="black" , xlim , ylim , buffer=c(0,0.4,0,0) , fix="" , skadoosh=TRUE , accel=-0.7 , vel0=1 , force_on="" , ... ) {

    

}
