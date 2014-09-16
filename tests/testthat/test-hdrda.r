library(testthat)
library(sparsediscrim)

context("The HDRDA Classifier from Ramey et al. (2014)")

test_that("The HDRDA classifier works properly on the iris data set", {
  require('MASS')

  set.seed(42)
  n <- nrow(iris)
  train <- sample(seq_len(n), n / 2)
  hdrda_out <- hdrda(Species ~ ., data = iris[train, ])
  predicted <- predict(hdrda_out, iris[-train, -5])$class

  hdrda_out2 <- hdrda(x = iris[train, -5], y = iris[train, 5])
  predicted2 <- predict(hdrda_out2, iris[-train, -5])$class

  # Tests that the same labels result from the matrix and formula versions of
  # the HDRDA classifier
  expect_equal(predicted, predicted2)
})

test_that("HDRDA's calculations are correct for (lambda, gamma) = (1, 0)", {
  require('MASS')

  hdrda_out <- hdrda(Species ~ ., data = iris, lambda=1, gamma=0)

  Sigma <- cov_pool(x=iris[, -5], y=iris$Species)
  Sigma_eigen <- eigen(Sigma, symmetric=TRUE)
  D_q <- diag(Sigma_eigen$values)

  # For each class, calculate W_k. These should equal D_q
  W_k <- sapply(hdrda_out$est, function(est) {
    solve(est$W_inv)
  }, simplify=FALSE)

  expect_equal(W_k$setosa, D_q)
  expect_equal(W_k$versicolor, D_q)
  expect_equal(W_k$virginica, D_q)
})

test_that("HDRDA's calculations are correct for (lambda, gamma) = (0, 0)", {
  require('MASS')

  hdrda_out <- hdrda(Species ~ ., data = iris, lambda=0, gamma=0)

  Sigma <- cov_pool(x=iris[, -5], y=iris$Species)
  Sigma_k <- tapply(seq_len(nrow(iris)), iris$Species, function(i) {
    cov_mle(iris[i, -5])
  })
  Sigma_eigen <- eigen(Sigma, symmetric=TRUE)
  U_1 <- Sigma_eigen$vectors

  # For each class, calculate W_k.
  W_k <- sapply(hdrda_out$est, function(est) {
    solve(est$W_inv)
  }, simplify=FALSE)

  W_k_expected <- sapply(Sigma_k, function(class_cov) {
    crossprod(U_1, class_cov) %*% U_1
  }, simplify=FALSE)

  expect_equal(W_k$setosa, W_k_expected$setosa)
  expect_equal(W_k$versicolor, W_k_expected$versicolor)
  expect_equal(W_k$virginica, W_k_expected$virginica)
})


test_that("LDA is a special case of HDRDA on the iris data set", {
  require('MASS')

  set.seed(42)
  n <- nrow(iris)
  train <- sample(seq_len(n), n / 2)
  hdrda_out <- hdrda(Species ~ ., data = iris[train, ], lambda=1, gamma=0)
  predicted_hdrda <- predict(hdrda_out, iris[-train, -5])$class

  lda_out <- lda(x = iris[train, -5], grouping = iris[train, 5])
  predicted_lda <- predict(lda_out, iris[-train, -5])$class

  # Tests that the same labels result from the matrix and formula versions of
  # the HDRDA classifier
  expect_equal(predicted_hdrda, predicted_lda)
})

test_that("QDA is a special case of HDRDA on the iris data set", {
  require('MASS')

  set.seed(42)
  n <- nrow(iris)
  train <- sample(seq_len(n), n / 2)
  hdrda_out <- hdrda(Species ~ ., data = iris[train, ], lambda=0, gamma=0)
  predicted_hdrda <- predict(hdrda_out, iris[-train, -5])$class

  qda_out <- qda(x = iris[train, -5], grouping = iris[train, 5])
  predicted_qda <- predict(qda_out, iris[-train, -5])$class

  # Tests that the same labels result from the matrix and formula versions of
  # the HDRDA classifier
  expect_equal(predicted_hdrda, predicted_qda)
})
