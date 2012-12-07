library('diagdiscrim')
library('testthat')

context("Simultaneous Diagonalization of Two Classes")

test_that("Two Low-Dimensional Covariance Matrices are Simultaneously Diagonalized with Rank Estimation to Select 'q'", {
	set.seed(42)
  # Sample size for each class
  n_k <- 20
	p <- 20
	x1 <- replicate(p, rnorm(n_k))
	x2 <- replicate(p, rnorm(n_k, 2))
  x <- rbind(x1, x2)
  y <- gl(2, n_k)

  simdiag_out <- simdiag(x, y)

  # Checks that the returned transformed matrix 'x' is the same as manually
  # transforming the data.
  transf_x1 <- tcrossprod(x1, simdiag_out$Q)
  transf_x2 <- tcrossprod(x2, simdiag_out$Q)
  transf_x <- rbind(transf_x1, transf_x2)
  expect_equal(simdiag_out$x, transf_x)

  # Checks that the transformed covariance matrices diagonal; numerically, they
  # are near diagonal -- the off-diagonal elements are "close" to zero but not
  # quite zero.
  cov_1 <- (n_k - 1) * cov(transf_x1) / n_k
  cov_2 <- (n_k - 1) * cov(transf_x2) / n_k
  expect_equal(cov_1, diag(diag(cov_1)))
  expect_equal(cov_2, diag(diag(cov_2)))
})

test_that("Two High-Dimensional Covariance Matrices are Simultaneously Diagonalized with Rank Estimation to Select 'q'", {
	set.seed(42)
  # Sample size for each class
  n_k <- 20
	p <- 100
	x1 <- replicate(p, rnorm(n_k))
	x2 <- replicate(p, rnorm(n_k, 2))
  x <- rbind(x1, x2)
  y <- gl(2, n_k)

  simdiag_out <- simdiag(x, y)

  # Checks that the returned transformed matrix 'x' is the same as manually
  # transforming the data.
  transf_x1 <- tcrossprod(x1, simdiag_out$Q)
  transf_x2 <- tcrossprod(x2, simdiag_out$Q)
  transf_x <- rbind(transf_x1, transf_x2)
  expect_equal(simdiag_out$x, transf_x)

  # Checks that the transformed covariance matrices diagonal; numerically, they
  # are near diagonal -- the off-diagonal elements are "close" to zero but not
  # quite zero.
  cov_1 <- (n_k - 1) * cov(transf_x1) / n_k
  cov_2 <- (n_k - 1) * cov(transf_x2) / n_k
  expect_equal(cov_1, diag(diag(cov_1)))
  expect_equal(cov_2, diag(diag(cov_2)))
})


