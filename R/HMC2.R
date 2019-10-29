# SIMPLE IMPLEMENTATION OF HAMILTONIAN MONTE CARLO.
# 
# Modified slightly based on code by 
# Radford M. Neal, 2010.
#
# This original appears in Figure 2 of "MCMC using Hamiltonian dynamics",
# the Handbook of Markov Chain Monte Carlo.
#
# The arguments to the HMC function are as follows:
#
#   U          A function to evaluate minus the log of the density of the
#              distribution to be sampled, plus any constant - ie, the
#              "potential energy".
#
#   grad_U     A function to evaluate the gradient of U.
#
#   epsilon    The stepsize to use for the leapfrog steps.
#
#   L          The number of leapfrog steps to do to propose a new state.
#
#   current_q  The current state (position variables only).
#
# Momentum variables are sampled from independent standard normal
# distributions within this function.  The value return is the vector
# of new position variables (equal to current_q if the endpoint of the
# trajectory was rejected).

# this modified version returns trajectories for q,p
# ... for arguments to pass to U and grad_U
# returns H for monitoring divergences
HMC2 <- function (U, grad_U, epsilon, L, current_q , ... ) {
  q = current_q
  p = rnorm(length(q),0,1)  # independent standard normal variates
  current_p = p
  # Make a half step for momentum at the beginning
  p = p - epsilon * grad_U( q , ... ) / 2
  # Alternate full steps for position and momentum
  qtraj <- matrix(NA,nrow=L+1,ncol=length(q))
  ptraj <- qtraj
  qtraj[1,] <- current_q
  ptraj[1,] <- p
  for (i in 1:L)
  {
    # Make a full step for the position
    q = q + epsilon * p
    # Make a full step for the momentum, except at end of trajectory
    if (i!=L) {
        p = p - epsilon * grad_U(q,...)
        ptraj[i+1,] <- p
    }
    qtraj[i+1,] <- q
  }
  # Make a half step for momentum at the end.
  p = p - epsilon * grad_U(q,...) / 2
  ptraj[L+1,] <- p
  # Negate momentum at end of trajectory to make the proposal symmetric
  p = -p
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = U(current_q,...)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q,...)
  proposed_K = sum(p^2) / 2
  H0 <- current_U + current_K
  H1 <- proposed_U + proposed_K
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  new_q <- q
  accept <- 0
  if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K)) {
    new_q <- q  # accept
    accept <- 1
  } else
    new_q <- current_q  # reject

  return(
    list(
        q = new_q,
        traj = qtraj,
        ptraj = ptraj,
        accept = accept , 
        dH = H1 - H0 ) )
}

HMC_2D_sample <- function( n=100 , U , U_gradient , step , L , start=c(0,0) , xlim=c(-5,5) , ylim=c(-4,4) , xlab="x" , ylab="y" , draw=TRUE , draw_contour=TRUE , nlvls=15 , adj_lvls=1 , ... ) {
    
    Q <- list()
    Q$q <- start
    xr <- xlim
    yr <- ylim

    if ( draw==TRUE ) plot( NULL , xlab=xlab , ylab=ylab , xlim=xlim, ylim=ylim )

    if ( draw==TRUE & draw_contour==TRUE ) {
        # draw contour
        zr <- 1
        y_seq <- seq(from=yr[1]-zr,to=yr[2]+zr,length.out=50) 
        x_seq <- seq(from=xr[1]-zr,to=xr[2]+zr,length.out=50) 
        z5 <- matrix(NA,length(x_seq),length(y_seq))
        for ( i in 1:length(x_seq) )
            for ( j in 1:length(y_seq) )
                z5[i,j] <- U( c( x_seq[i] , y_seq[j] ) , ... )
        lz5 <- log(z5)
        lvls <- exp( pretty( range(lz5), nlvls ) )
        cl <- contourLines( x_seq , y_seq , z5 , level=lvls*adj_lvls )
        for ( i in 1:length(cl) ) lines( cl[[i]]$x , cl[[i]]$y , col=col.alpha("black",0.6) , lwd=0.5 )
    }

    n_samples <- n
    a <- rep(NA,n_samples)
    dH <- rep(NA,n_samples) # energy
    path_col <- col.alpha("black",0.3)
    xpos <- c(1,3,4,2)
    points( Q$q[1] , Q$q[2] , pch=4 , col="black" )
    post <- matrix(NA,nrow=n,ncol=2)

    for ( i in 1:n_samples ) {
        Q <- HMC2( U , U_gradient , epsilon=step , L=L , current_q=Q$q )
        if ( n_samples < 100 ) {
          # draw paths
          lines( Q$traj[,1] , Q$traj[,2] , col=path_col , lwd=2 )
        }
        r <- min(abs(Q$dH),1)
        ptcol <- rgb( r , 0 , 0 )
        ptcol <- "black"
        points( Q$traj[L+1,1] , Q$traj[L+1,2] , pch=ifelse( Q$accept==1 , 16 , 1 ) , col=ptcol , cex=ifelse( Q$accept==1 , 0.7 , 1 ) , lwd=1.5 )
        dH[i] <- Q$dH
        a[i] <- Q$accept
        if ( a[i]==1 ) {
            post[i,] <- Q$q
        }
    }
    
    return( invisible( post ) )
}
