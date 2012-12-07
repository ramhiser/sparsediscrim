library('diagdiscrim')
library('testthat')

library('plsgenomics')

context("Simultaneous Diagonalization of the Alon Data Set with the Fast SVD")

data('Colon')
x <- Colon$X
y <- Colon$Y

y <- factor(y)
levels_y <- levels(y)

x1 <- x[which(y == levels_y[1]), ]
x2 <- x[which(y == levels_y[2]), ]

cov_x1 <- diagdiscrim:::cov_mle(x1)
cov_x2 <- diagdiscrim:::cov_mle(x2)

# By default, the first class is chosen to be the one with the larger sample size.
# The reason is that if there is no multicollinearity, then the class with larger
# sample size should have a higher rank sample covariance matrix (assuming that
# the data are small versus the dimension of the data). Our formulation is set up
# so that the first class should be  the one that has a higher rank sample
# covariance matrix.
test_that("The Alon Data Set is Simultaneously Diagonalized based on Sample Size", {
  Q <- diagdiscrim:::simdiag_fastsvd(x = x, y = y)$Q

  simdiag_cov_1 <- Q %*% cov_x1 %*% t(Q)
  simdiag_cov_2 <- Q %*% cov_x2 %*% t(Q)  
  expect_equal(simdiag_cov_1, diag(diag(simdiag_cov_1)))
  expect_equal(simdiag_cov_2, diag(diag(simdiag_cov_2)))
})

# If multicollinearity is suspected and if the sample sizes for both classes are
# approximately equal, then it may not be the case that the class with larger
# sample size should have a higher rank sample covariance matrix. In this case,
# we provide the user with the option to numerically estimate the ranks of each
# class and choose the first class to be the one with a larger estimated rank.
# This is can be more more accurate than the default tested above, but also more
# computationally intense because the "Fast SVD" is performed on both classes.
test_that("The Alon Data Set is Simultaneously Diagonalized based on Estimated Larger Rank", {
  Q <- diagdiscrim:::simdiag_fastsvd(x = x, y = y, calc_ranks = TRUE)$Q

  simdiag_cov_1 <- Q %*% cov_x1 %*% t(Q)
  simdiag_cov_2 <- Q %*% cov_x2 %*% t(Q)  
  expect_equal(simdiag_cov_1, diag(diag(simdiag_cov_1)))
  expect_equal(simdiag_cov_2, diag(diag(simdiag_cov_2)))
})

