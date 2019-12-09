\name{HMC2}
\alias{HMC2}\alias{HMC_2D_sample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Functions for simple HMC simulations}
\description{
  Conduct simple Hamiltonian Monte Carlo simulations.
}
\usage{
HMC2(U, grad_U, epsilon, L, current_q , ... )
HMC_2D_sample( n=100 , U , U_gradient , step , L , start=c(0,0) , 
  xlim=c(-5,5) , ylim=c(-4,4) , xlab="x" , ylab="y" , 
  draw=TRUE , draw_contour=TRUE , nlvls=15 , adj_lvls=1 , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{U}{Function to return log density}
  \item{grad_U}{Function to return gradient of U}
  \item{epsilon}{Size of leapfrog step}
  \item{L}{Number of leapfrog steps}
  \item{current_q}{Initial position}
  \item{n}{Number of transitions to sample}
  \item{step}{Size of leapfrog step}
  \item{start}{Initial position}
  \item{xlim}{Horizontal boundaries of plot, for when \code{draw} is TRUE}
  \item{ylim}{Vertical boundaries of plot, for when \code{draw} is TRUE}
  \item{draw}{Whether to draw the samples and trajectories}
  \item{draw_contour}{Whether to draw contour of density}
  \item{nlvls}{Number of contour levels}
  \item{adj_lvls}{Factor to multiply levels by}
  \item{...}{Optional arguments to pass to density and gradient functions, typically optional parameters}
}
\details{
  These functions provide simple Hamiltonian Monte Carlo simulations.

  \code{HMC2} is based on Neals's 2010 code (see citation below), but returns the trajectories.

  \code{HMC_2D_sample} simulates multiple sequential trajectories and can also plot the simulation. See examples below.

  To use either function, the user must supply custom density and gradient functions.
}
\value{
}
\references{}
\author{Richard McElreath}
\seealso{Radford M. Neal, 2010. MCMC using Hamiltonian dynamics. The Handbook of Markov Chain Monte Carlo.}
\examples{
# Devil's Funnel from Chapter 13
U_funnel <- function( q , s=3 ) {
    v <- q[2]
    x <- q[1]
    U <- sum( dnorm(x,0,exp(v),log=TRUE) ) + dnorm(v,0,s,log=TRUE)
    return( -U )
}
U_funnel_gradient <- function( q , s=3 ) {
    v <- q[2]
    x <- q[1]
    Gv <- (-v)/s^2 - length(x) + exp(-2*v)*sum( x^2 ) #dU/dv
    Gx <- -exp(-2*v)*x #dU/dx
    return( c( -Gx , -Gv ) ) # negative bc energy is neg-log-prob
}
HMC_2D_sample( n=3 , U=U_funnel , U_gradient=U_funnel_gradient , 
  step=0.2 , L=10 , ylab="v"  , adj_lvls=1/12 )

# Same, but with non-centered parameterization

U_funnel_NC <- function( q , s=3 ) {
    v <- q[2]
    z <- q[1]
    U <- sum( dnorm(z,0,1,log=TRUE) ) + dnorm(v,0,s,log=TRUE)
    return( -U )
}
U_funnel_NC_gradient <- function( q , s=3 ) {
    v <- q[2]
    z <- q[1]
    Gv <- (-v)/s^2 #dU/dv
    Gz <- (-z) #dU/dz
    return( c( -Gz , -Gv ) ) # negative bc energy is neg-log-prob
}
HMC_2D_sample( n=4 , U=U_funnel_NC , U_gradient=U_funnel_NC_gradient , 
  step=0.2 , L=15 , ylab="v" , xlab="z" , xlim=c(-2,2) , nlvls=12 , adj_lvls=1 )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

