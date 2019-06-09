\name{trankplot}
\alias{trankplot}
\alias{traceplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Diagnostic trace and rank histogram plots for MCMC output}
\description{
  The functions \code{trankplot} and \code{traceplot} display MCMC chain diagnostic plots. \code{trankplot} displays ranked histograms and \code{traceplot} shows the more traditional trace of the samples.
}
\usage{
trankplot( object , bins=30 , pars , chains , col=rethink_palette , alpha=1 , 
  bg=col.alpha("black",0.15) , ask=TRUE , window , n_cols=3 , max_rows=5 , 
  lwd=1.5 , lp=FALSE  , axes=FALSE , ... )
traceplot( object , pars , chains , col=rethink_palette , alpha=1 , 
  bg=col.alpha("black",0.15) , ask=TRUE , window , trim=100 , n_cols=3 , 
  max_rows=5 , lwd=0.5 , lp=FALSE , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A \code{stanfit}, \code{ulam} or \code{map2stan} object}
  \item{bins}{For \code{trankplot}, the number of histogram bins to use}
  \item{pars}{Optional character vector of parameters to display}
  \item{chains}{Optional integer vector of chains to display}
  \item{col}{Vector of colors to use for chains}
  \item{alpha}{Transparency}
  \item{bg}{Background color for warmup samples}
  \item{ask}{Interactive paging when \code{TRUE}}
  \item{window}{Optional range of samples to show}
  \item{n_cols}{Number of columns in display}
  \item{max_rows}{Maximum number of rows on each page}
  \item{lwd}{Line width}
  \item{lp}{Whether to include log_prob in display}
  \item{axes}{Whether to show axes on plots}
  \item{trim}{For \code{traceplot}, number of samples to trim for start. Helps with display, since early warmup samples typically very far from typical samples.}
  \item{...}{Additional arguments to pass to \code{\link{plot}}}
}
\details{
  \code{trankplot} produces rank histograms of each chain, as described in Vehtari et al 2019 (see reference below). For each parameter, the samples from all chains are first ranked, using \code{rank_mat}. This returns a matrix of ranks, with the chains preserved. Then a histogram is built for each chain, using the same break points. These historgrams are then overlain in the plot.

  For healthy well-mixing chains, the histrograms should be uniform. When there are spikes for some chains, especially in the low or high ranks, this suggests problems in exploring the posterior.

  \code{traceplot} shows the sequential samples for each parameter and chain. This is the same information as the \code{trankplot}, but often much harder to see, given the volume of samples.
}
\references{Vehtari, Gelman, Simpson, Carpenter, Bürkner. 2019. Rank-normalization, folding, and localization: An improved R-hat for assessing convergence of MCMC. https://arxiv.org/abs/1903.08008}
\author{Richard McElreath}
\seealso{}
\examples{
\dontrun{
library(rethinking)
data(chimpanzees)

d <- list( 
    pulled_left = chimpanzees$pulled_left ,
    prosoc_left = chimpanzees$prosoc_left ,
    condition = as.integer( 2 - chimpanzees$condition ) ,
    actor = as.integer( chimpanzees$actor ) ,
    blockid = as.integer( chimpanzees$block )
)

m <- ulam(
    alist(
        # likeliood
        pulled_left ~ bernoulli(theta),

        # linear models
        logit(theta) <- A + BP*prosoc_left,
        A <- a + v[actor,1],
        BP <- bp + v[actor,condition+1],

        # adaptive prior
        vector[3]: v[actor] ~ multi_normal( 0 , Rho_actor , sigma_actor ),

        # fixed priors
        c(a,bp) ~ normal(0,1),
        sigma_actor ~ exponential(1),
        Rho_actor ~ lkjcorr(4)
    ) , data=d , chains=3 , cores=1 , sample=TRUE )

trankplot(m)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

