\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}

\name{NEWS}
\title{rethinking News}
\encoding{UTF-8}
\section{CHANGES IN VERSION 1.2.7 (2013-12-29)}{
  \itemize{
     \item added \code{map} convergence warnings
     \item added recognition of log and logit links to \code{map} 
     \item \code{sample.qa.posterior} can now accept arbitrary model weights
   }
}
\section{CHANGES IN VERSION 1.2.0 (2013-07-11)}{
  \itemize{
     \item add \code{map} command to replace functionality of \code{mle2} in \code{bbmle} package
     \item added NEWS (this document)
   }
}