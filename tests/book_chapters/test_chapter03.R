# test_dir('rethinking/tests/book_chapters',reporter="summary")
# test_dir('rethinking/tests/book_chapters',filter="chapter03")

context('chapter 3')
library(rethinking)

expect_equiv_eps <- function( x , y , eps=0.01 ) {
    expect_equivalent( x , y , tolerance=eps )
}

set.seed(39487)

## R code 3.11
p_grid <- seq( from=0 , to=1 , length.out=1000 )
prior <- rep(1,1000)
likelihood <- dbinom( 3 , size=3 , prob=p_grid )
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)
samples <- sample( p_grid , size=1e4 , replace=TRUE , prob=posterior )

## R code 3.12
xxx <- PI( samples , prob=0.5 )

test_that("R code 3.12",
    expect_equiv_eps( round(xxx,1) , c(0.7,0.9) )
)

## R code 3.13
rm(xxx)
xxx <- HPDI( samples , prob=0.5 )

test_that("R code 3.13",
    expect_equiv_eps( round(xxx,2) , c(0.84 , 1 ) )
)
