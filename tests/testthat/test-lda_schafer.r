library(testthat)
library(sparsediscrim)

context("LDA with the Schafer-Strimmer Covariance Matrix")

test_that("The lda_schafer classifier works properly on the iris data set", {
  require('MASS')

  set.seed(42)
  n <- nrow(iris)
  train <- sample(seq_len(n), n / 2)
  lda_schafer_out <- lda_schafer(Species ~ ., data = iris[train, ])
  predicted <- predict(lda_schafer_out, iris[-train, -5])

  lda_schafer_out2 <- lda_schafer(x = iris[train, -5], y = iris[train, 5])
  predicted2 <- predict(lda_schafer_out2, iris[-train, -5])

  # Tests that the same labels result from the matrix and formula versions of
  # the lda_schafer classifier
  expect_equal(predicted$class, predicted2$class)

  expect_is(predicted$posterior, "matrix")
})

# Related to issue #41
test_that("The lda_schafer classifier works properly when 1 feature used", {
  require('MASS')

  set.seed(42)
  n <- nrow(iris)
  train <- sample(seq_len(n), n / 2)
  n_test <- n - length(train)

  lda_schafer_out <- lda_schafer(x = iris[train, 1], y = iris[train, 5])
  predicted <- predict(lda_schafer_out, iris[-train, 1])

  expect_equal(length(predicted$class), n_test)
  expect_is(predicted$posterior, "matrix")
  expect_equal(dim(predicted$posterior), c(n_test, nlevels(iris$Species)))
})
