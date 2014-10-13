library(testthat)
library(sparsediscrim)

context("Generate multivariate normal variates with block-diagonal autocorrelation matrices")

test_that("Autocorrelation matrix is generated properly", {
  cov_matrix <- cov_autocorrelation(p=3, rho=0.9)
  cov_expected <- cbind(c(1, 0.9, 0.81), c(0.9, 1, 0.9), c(0.81, 0.9, 1))

  expect_equal(cov_matrix, cov_expected)
})

test_that("Block-diagonal covariance matrix with autocorrelated blocks is generated properly", {
  cov_matrix <- cov_block_autocorrelation(num_blocks=2, block_size=3, rho=0.9)

  autocorr_block <- cbind(c(1, 0.9, 0.81), c(0.9, 1, 0.9), c(0.81, 0.9, 1))
  zeros_matrix <- matrix(0, nrow=3, ncol=3)

  cov_expected1 <- cbind(autocorr_block, zeros_matrix)
  cov_expected2 <- cbind(zeros_matrix, autocorr_block)
  cov_expected <- rbind(cov_expected1, cov_expected2)

  expect_equal(cov_matrix, cov_expected, check.attributes=FALSE)
})

test_that("Generate MVN data set with block-diagonal autocorrelation matrices", {
  n <- c(10, 10, 10)
  means <- matrix(rep(1:3, each=6), ncol=3)
  cov_matrix <- cov_block_autocorrelation(num_blocks=2, block_size=3, rho=0.9)

  # Generated Data
  set.seed(42)
  data_gen <- generate_blockdiag(n=n, mu=means, num_blocks=2, block_size=3, rho=rep(0.9, 3))

  # Expected Data
  set.seed(42)
  x1 <- rmvnorm(n=n[1], mean=means[, 1], sigma=cov_matrix)
  x2 <- rmvnorm(n=n[2], mean=means[, 2], sigma=cov_matrix)
  x3 <- rmvnorm(n=n[3], mean=means[, 3], sigma=cov_matrix)
  x_expected <- rbind(x1, x2, x3)
  y_expected <- factor(rep(seq_along(n), n))

  expect_equal(data_gen$x, x_expected)
  expect_equal(data_gen$y, y_expected)
})
