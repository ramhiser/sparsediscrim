library(testthat)
library(sparsediscrim)

context("The DLDA Classifier from Dudoit et al. (2002)")

test_that("The DLDA classifier works properly on the iris data set", {
  require('MASS')

  set.seed(42)
  n <- nrow(iris)
  train <- sample(seq_len(n), n / 2)
  dlda_out <- lda_diag(Species ~ ., data = iris[train, ])
  predicted <- predict(dlda_out, iris[-train, -5])

  dlda_out2 <- lda_diag(x = iris[train, -5], y = iris[train, 5])
  predicted2 <- predict(dlda_out2, iris[-train, -5])

  # Tests that the same labels result from the matrix and formula versions of
  # the DLDA classifier
  expect_equal(predicted$class, predicted2$class)

  expect_is(predicted$posterior, "matrix")
})

# Related to issue #41
test_that("The DLDA classifier works properly when 1 feature used", {
  require('MASS')

  set.seed(42)
  n <- nrow(iris)
  train <- sample(seq_len(n), n / 2)
  n_test <- n - length(train)

  dlda_out <- lda_diag(x = iris[train, 1], y = iris[train, 5])
  predicted <- predict(dlda_out, iris[-train, 1])

  expect_equal(length(predicted$class), n_test)
  expect_is(predicted$posterior, "matrix")
  expect_equal(dim(predicted$posterior), c(n_test, nlevels(iris$Species)))
})
