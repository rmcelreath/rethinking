# test_dir('rethinking/tests/ulam_tests',reporter="summary")

context('ulam')
library(rethinking)
# library(openssl) # for md5 hash

expect_equiv_eps <- function( x , y , eps=0.01 ) {
    expect_equivalent( x , y , tolerance=eps )
}

# student_t

data(Howell1)
d <- Howell1
d2 <- d[ d$age >= 18 , ]

set.seed(1)
mt <- ulam(
  alist(
    height ~ student_t( nu, mu , sigma ) ,
    mu <- a + b*weight ,
    a ~ normal( 178 , 20 ) ,
    b ~ normal( 0 , 1 ) ,
    sigma ~ exponential( 1 ),
    nu ~ gamma(2, 0.1)
  ) , data=d2 , sample=TRUE , refresh=-1 , seed=1 )

test_that("Student_t regression code match",
    expect_known_hash( mt@model , "98532fc851" )
)

coef_expect <- c(113.62  , 0.91  , 3.96  , 3.42  )
names(coef_expect) <- c("a","b","sigma","nu")
test_that("Student_t regression post means",
    expect_equiv_eps( round(coef(mt),2) , coef_expect )
)

rm(mt)
rm(d)
rm(d2)
gc()

# pareto

set.seed(49283)
y <- rpareto(1e3,1,1)
mp <- ulam(
    alist(
        y ~ pareto(1,alpha),
        alpha ~ exponential(1)
    ) , data=list(y=y) , chains=1 )

test_that("Pareto regression code match",
    expect_known_hash( mp@model , "7d9808c554" )
)

test_that("Pareto regression post match",
    expect_equiv_eps( round(coef(mp),2) , 1.01 )
)


