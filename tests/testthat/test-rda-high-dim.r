library(testthat)
library(sparsediscrim)
library(MASS)
library(mvtnorm)

data(two_class_sim_data)

context("The HDRDA Classifier from Ramey et al. (2017)")

test_that("The HDRDA classifier works properly on the iris data set", {
  set.seed(42)
  n <- nrow(iris)
  train <- sample(seq_len(n), n / 2)
  hdrda_out <- rda_high_dim(Species ~ ., data = iris[train, ])
  predicted <- predict(hdrda_out, iris[-train, -5])$class

  hdrda_out2 <- rda_high_dim(x = iris[train, -5], y = iris[train, 5])
  predicted2 <- predict(hdrda_out2, iris[-train, -5])$class

  # Tests that the same labels result from the matrix and formula versions of
  # the HDRDA classifier
  expect_equal(predicted, predicted2)
})

test_that("HDRDA's calculations are correct for (lambda, gamma) = (1, 0)", {
  hdrda_out <- rda_high_dim(Species ~ ., data = iris, lambda=1, gamma=0)

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

test_that("LDA is a special case of HDRDA on the iris data set", {
  set.seed(42)
  n <- nrow(iris)
  train <- sample(seq_len(n), n / 2)
  hdrda_out <- rda_high_dim(Species ~ ., data = iris[train, ], lambda=1, gamma=0)
  predicted_hdrda <- predict(hdrda_out, iris[-train, -5])$class

  lda_out <- lda(x = iris[train, -5], grouping = iris[train, 5])
  predicted_lda <- predict(lda_out, iris[-train, -5])$class

  # Tests that the same labels result from the matrix and formula versions of
  # the HDRDA classifier
  expect_equal(predicted_hdrda, predicted_lda)
})

test_that("QDA is a special case of HDRDA on the iris data set", {
  set.seed(42)
  n <- nrow(iris)
  train <- sample(seq_len(n), n / 2)
  hdrda_out <- rda_high_dim(Species ~ ., data = iris[train, ], lambda=0, gamma=0)
  predicted_hdrda <- predict(hdrda_out, iris[-train, -5])$class

  qda_out <- qda(x = iris[train, -5], grouping = iris[train, 5])
  predicted_qda <- predict(qda_out, iris[-train, -5])$class

  # Tests that the same labels result from the matrix and formula versions of
  # the HDRDA classifier
  expect_equal(predicted_hdrda, predicted_qda)
})

test_that("HDRDA's statistics match manual values when (lambda, gamma) = (0, 0)", {
  set.seed(42)
  p <- 250
  n1 <- 10
  n2 <- 15
  n3 <- 25

  iris_x <- as.matrix(iris[, -5])

  # The first cov matrix is the identity matrix.
  # To generate sigma2 and sigma3, we grab eigenvectors from the iris data set.
  # The eigenvalues are randomly generated.
  eigen2 <- eigen(tcrossprod(iris_x[rep(1:100, length=p), ]), symmetric=TRUE)
  eigen3 <- eigen(tcrossprod(iris_x[rep(30:130, length=p), ]), symmetric=TRUE)

  eigen2$values <- sort(rexp(p, rate=0.1), decreasing=TRUE)
  eigen3$values <- sort(rexp(p, rate=0.05), decreasing=TRUE)

  sigma2 <- with(eigen2, vectors %*% tcrossprod(diag(values), vectors))
  sigma3 <- with(eigen3, vectors %*% tcrossprod(diag(values), vectors))

  x1 <- rmvnorm(n1, mean=rep(0, p), sigma=diag(p))
  x2 <- rmvnorm(n2, mean=rep(3, p), sigma=sigma2)
  x3 <- rmvnorm(n3, mean=rep(-3, p), sigma=sigma3)

  x <- rbind(x1, x2, x3)
  y <- c(rep(1, n1), rep(2, n2), rep(3, n3))

  # Compute estimates for HDRDA manually.
  xbar1 <- colMeans(x1)
  xbar2 <- colMeans(x2)
  xbar3 <- colMeans(x3)

  x_centered1 <- scale(x1, center=TRUE, scale=FALSE)
  x_centered2 <- scale(x2, center=TRUE, scale=FALSE)
  x_centered3 <- scale(x3, center=TRUE, scale=FALSE)

  Sigma_eigen <- cov_eigen(x=x, y=y, pool=TRUE, fast=TRUE, tol=1e-6)
  D_q <- Sigma_eigen$values
  U1 <- Sigma_eigen$vectors
  q <- length(D_q)

  XU_1 <- x_centered1 %*% U1
  XU_2 <- x_centered2 %*% U1
  XU_3 <- x_centered3 %*% U1

  xbar1_U1 <- crossprod(U1, xbar1)
  xbar2_U1 <- crossprod(U1, xbar2)
  xbar3_U1 <- crossprod(U1, xbar3)

  W_1 <- cov_mle(XU_1)
  W_2 <- cov_mle(XU_2)
  W_3 <- cov_mle(XU_3)

  W1_inv <- solve_chol(W_1 + diag(0.001, nrow=nrow(W_1), ncol=ncol(W_1)))
  W2_inv <- solve_chol(W_2 + diag(0.001, nrow=nrow(W_2), ncol=ncol(W_2)))
  W3_inv <- solve_chol(W_3 + diag(0.001, nrow=nrow(W_3), ncol=ncol(W_3)))

  hdrda_out <- rda_high_dim(x=x, y=y, lambda=0, gamma=0)

  expect_equal(hdrda_out$q, q)
  expect_equal(hdrda_out$D_q, D_q)
  expect_equal(hdrda_out$U1, U1)

  expect_equal(hdrda_out$est[[1]]$n, n1)
  expect_equal(hdrda_out$est[[1]]$xbar, xbar1)
  expect_equal(hdrda_out$est[[1]]$prior, n1 / length(y))
  expect_equal(hdrda_out$est[[1]]$alpha, 1)
  expect_equal(hdrda_out$est[[1]]$XU, XU_1)
  expect_equal(hdrda_out$est[[1]]$xbar_U1, xbar1_U1)
  expect_equal(hdrda_out$est[[1]]$Gamma, matrix(0, q, q))
  expect_equal(hdrda_out$est[[1]]$Q, diag(n1))
  expect_equal(hdrda_out$est[[1]]$W_inv, W1_inv)

  expect_equal(hdrda_out$est[[2]]$n, n2)
  expect_equal(hdrda_out$est[[2]]$xbar, xbar2)
  expect_equal(hdrda_out$est[[2]]$prior, n2 / length(y))
  expect_equal(hdrda_out$est[[2]]$alpha, 1)
  expect_equal(hdrda_out$est[[2]]$XU, XU_2)
  expect_equal(hdrda_out$est[[2]]$xbar_U1, xbar2_U1)
  expect_equal(hdrda_out$est[[2]]$Gamma, matrix(0, q, q))
  expect_equal(hdrda_out$est[[2]]$Q, diag(n2))
  expect_equal(hdrda_out$est[[2]]$W_inv, W2_inv)

  expect_equal(hdrda_out$est[[3]]$n, n3)
  expect_equal(hdrda_out$est[[3]]$xbar, xbar3)
  expect_equal(hdrda_out$est[[3]]$prior, n3 / length(y))
  expect_equal(hdrda_out$est[[3]]$alpha, 1)
  expect_equal(hdrda_out$est[[3]]$XU, XU_3)
  expect_equal(hdrda_out$est[[3]]$xbar_U1, xbar3_U1)
  expect_equal(hdrda_out$est[[3]]$Gamma, matrix(0, q, q))
  expect_equal(hdrda_out$est[[3]]$Q, diag(n3))
  expect_equal(hdrda_out$est[[3]]$W_inv, W3_inv)
})

test_that("HDRDA's statistics match manual values when (lambda, gamma) = (1, 0)", {
  set.seed(42)
  p <- 250
  n1 <- 10
  n2 <- 15
  n3 <- 25

  iris_x <- as.matrix(iris[, -5])

  # The first cov matrix is the identity matrix.
  # To generate sigma2 and sigma3, we grab eigenvectors from the iris data set.
  # The eigenvalues are randomly generated.
  eigen2 <- eigen(tcrossprod(iris_x[rep(1:100, length=p), ]), symmetric=TRUE)
  eigen3 <- eigen(tcrossprod(iris_x[rep(30:130, length=p), ]), symmetric=TRUE)

  eigen2$values <- sort(rexp(p, rate=0.1), decreasing=TRUE)
  eigen3$values <- sort(rexp(p, rate=0.05), decreasing=TRUE)

  sigma2 <- with(eigen2, vectors %*% tcrossprod(diag(values), vectors))
  sigma3 <- with(eigen3, vectors %*% tcrossprod(diag(values), vectors))

  x1 <- rmvnorm(n1, mean=rep(0, p), sigma=diag(p))
  x2 <- rmvnorm(n2, mean=rep(3, p), sigma=sigma2)
  x3 <- rmvnorm(n3, mean=rep(-3, p), sigma=sigma3)

  x <- rbind(x1, x2, x3)
  y <- c(rep(1, n1), rep(2, n2), rep(3, n3))

  # Compute estimates for HDRDA manually.
  xbar1 <- colMeans(x1)
  xbar2 <- colMeans(x2)
  xbar3 <- colMeans(x3)

  x_centered1 <- scale(x1, center=TRUE, scale=FALSE)
  x_centered2 <- scale(x2, center=TRUE, scale=FALSE)
  x_centered3 <- scale(x3, center=TRUE, scale=FALSE)

  Sigma_eigen <- cov_eigen(x=x, y=y, pool=TRUE, fast=TRUE, tol=1e-6)
  D_q <- Sigma_eigen$values
  U1 <- Sigma_eigen$vectors
  q <- length(D_q)

  XU_1 <- x_centered1 %*% U1
  XU_2 <- x_centered2 %*% U1
  XU_3 <- x_centered3 %*% U1

  xbar1_U1 <- crossprod(U1, xbar1)
  xbar2_U1 <- crossprod(U1, xbar2)
  xbar3_U1 <- crossprod(U1, xbar3)

  hdrda_out <- rda_high_dim(x=x, y=y, lambda=1, gamma=0)

  expect_equal(hdrda_out$q, q)
  expect_equal(hdrda_out$D_q, D_q)
  expect_equal(hdrda_out$U1, U1)

  expect_equal(hdrda_out$est[[1]]$n, n1)
  expect_equal(hdrda_out$est[[1]]$xbar, xbar1)
  expect_equal(hdrda_out$est[[1]]$prior, n1 / length(y))
  expect_equal(hdrda_out$est[[1]]$alpha, 1)
  expect_equal(hdrda_out$est[[1]]$XU, XU_1)
  expect_equal(hdrda_out$est[[1]]$xbar_U1, xbar1_U1)
  expect_equal(hdrda_out$est[[1]]$Gamma, D_q)
  expect_equal(hdrda_out$est[[1]]$Q, diag(n1))
  expect_equal(hdrda_out$est[[1]]$W_inv, diag(D_q^(-1)))

  expect_equal(hdrda_out$est[[2]]$n, n2)
  expect_equal(hdrda_out$est[[2]]$xbar, xbar2)
  expect_equal(hdrda_out$est[[2]]$prior, n2 / length(y))
  expect_equal(hdrda_out$est[[2]]$alpha, 1)
  expect_equal(hdrda_out$est[[2]]$XU, XU_2)
  expect_equal(hdrda_out$est[[2]]$xbar_U1, xbar2_U1)
  expect_equal(hdrda_out$est[[2]]$Gamma, D_q)
  expect_equal(hdrda_out$est[[2]]$Q, diag(n2))
  expect_equal(hdrda_out$est[[2]]$W_inv, diag(D_q^(-1)))

  expect_equal(hdrda_out$est[[3]]$n, n3)
  expect_equal(hdrda_out$est[[3]]$xbar, xbar3)
  expect_equal(hdrda_out$est[[3]]$prior, n3 / length(y))
  expect_equal(hdrda_out$est[[3]]$alpha, 1)
  expect_equal(hdrda_out$est[[3]]$XU, XU_3)
  expect_equal(hdrda_out$est[[3]]$xbar_U1, xbar3_U1)
  expect_equal(hdrda_out$est[[3]]$Gamma, D_q)
  expect_equal(hdrda_out$est[[3]]$Q, diag(n3))
  expect_equal(hdrda_out$est[[3]]$W_inv, diag(D_q^(-1)))
})

test_that("HDRDA's statistics match manual values when (lambda, gamma) = (0.5, 0.5)", {
  set.seed(42)
  p <- 250
  n1 <- 10
  n2 <- 15
  n3 <- 25

  iris_x <- as.matrix(iris[, -5])

  # The first cov matrix is the identity matrix.
  # To generate sigma2 and sigma3, we grab eigenvectors from the iris data set.
  # The eigenvalues are randomly generated.
  eigen2 <- eigen(tcrossprod(iris_x[rep(1:100, length=p), ]), symmetric=TRUE)
  eigen3 <- eigen(tcrossprod(iris_x[rep(30:130, length=p), ]), symmetric=TRUE)

  eigen2$values <- sort(rexp(p, rate=0.1), decreasing=TRUE)
  eigen3$values <- sort(rexp(p, rate=0.05), decreasing=TRUE)

  sigma2 <- with(eigen2, vectors %*% tcrossprod(diag(values), vectors))
  sigma3 <- with(eigen3, vectors %*% tcrossprod(diag(values), vectors))

  x1 <- rmvnorm(n1, mean=rep(0, p), sigma=diag(p))
  x2 <- rmvnorm(n2, mean=rep(3, p), sigma=sigma2)
  x3 <- rmvnorm(n3, mean=rep(-3, p), sigma=sigma3)

  x <- rbind(x1, x2, x3)
  y <- c(rep(1, n1), rep(2, n2), rep(3, n3))

  # Compute estimates for HDRDA manually.
  xbar1 <- colMeans(x1)
  xbar2 <- colMeans(x2)
  xbar3 <- colMeans(x3)

  x_centered1 <- scale(x1, center=TRUE, scale=FALSE)
  x_centered2 <- scale(x2, center=TRUE, scale=FALSE)
  x_centered3 <- scale(x3, center=TRUE, scale=FALSE)

  Sigma_eigen <- cov_eigen(x=x, y=y, pool=TRUE, fast=TRUE, tol=1e-6)
  D_q <- Sigma_eigen$values
  U1 <- Sigma_eigen$vectors
  q <- length(D_q)

  XU_1 <- x_centered1 %*% U1
  XU_2 <- x_centered2 %*% U1
  XU_3 <- x_centered3 %*% U1

  xbar1_U1 <- crossprod(U1, xbar1)
  xbar2_U1 <- crossprod(U1, xbar2)
  xbar3_U1 <- crossprod(U1, xbar3)

  Gamma <- 0.5 * diag(D_q) + 0.5 * diag(q)
  Gamma_inv <- diag(diag(Gamma)^(-1))

  W_1 <- 0.5 * cov_mle(XU_1) + Gamma
  W_2 <- 0.5 * cov_mle(XU_2) + Gamma
  W_3 <- 0.5 * cov_mle(XU_3) + Gamma

  W1_inv <- solve_chol(W_1)
  W2_inv <- solve_chol(W_2)
  W3_inv <- solve_chol(W_3)

  Q1 <- diag(n1) + 0.5 / n1 * XU_1 %*% tcrossprod(Gamma_inv, XU_1)
  Q2 <- diag(n2) + 0.5 / n2 * XU_2 %*% tcrossprod(Gamma_inv, XU_2)
  Q3 <- diag(n3) + 0.5 / n3 * XU_3 %*% tcrossprod(Gamma_inv, XU_3)

  hdrda_out <- rda_high_dim(x=x, y=y, lambda=0.5, gamma=0.5)

  expect_equal(hdrda_out$q, q)
  expect_equal(hdrda_out$D_q, D_q)
  expect_equal(hdrda_out$U1, U1)

  expect_equal(hdrda_out$est[[1]]$n, n1)
  expect_equal(hdrda_out$est[[1]]$xbar, xbar1)
  expect_equal(hdrda_out$est[[1]]$prior, n1 / length(y))
  expect_equal(hdrda_out$est[[1]]$alpha, 1)
  expect_equal(hdrda_out$est[[1]]$XU, XU_1)
  expect_equal(hdrda_out$est[[1]]$xbar_U1, xbar1_U1)
  expect_equal(hdrda_out$est[[1]]$Gamma, diag(Gamma))
  expect_equal(hdrda_out$est[[1]]$Q, Q1)
  expect_equal(hdrda_out$est[[1]]$W_inv, W1_inv)

  expect_equal(hdrda_out$est[[2]]$n, n2)
  expect_equal(hdrda_out$est[[2]]$xbar, xbar2)
  expect_equal(hdrda_out$est[[2]]$prior, n2 / length(y))
  expect_equal(hdrda_out$est[[2]]$alpha, 1)
  expect_equal(hdrda_out$est[[2]]$XU, XU_2)
  expect_equal(hdrda_out$est[[2]]$xbar_U1, xbar2_U1)
  expect_equal(hdrda_out$est[[2]]$Gamma, diag(Gamma))
  expect_equal(hdrda_out$est[[2]]$Q, Q2)
  expect_equal(hdrda_out$est[[2]]$W_inv, W2_inv)

  expect_equal(hdrda_out$est[[3]]$n, n3)
  expect_equal(hdrda_out$est[[3]]$xbar, xbar3)
  expect_equal(hdrda_out$est[[3]]$prior, n3 / length(y))
  expect_equal(hdrda_out$est[[3]]$alpha, 1)
  expect_equal(hdrda_out$est[[3]]$XU, XU_3)
  expect_equal(hdrda_out$est[[3]]$xbar_U1, xbar3_U1)
  expect_equal(hdrda_out$est[[3]]$Gamma, diag(Gamma))
  expect_equal(hdrda_out$est[[3]]$Q, Q3)
  expect_equal(hdrda_out$est[[3]]$W_inv, W3_inv)
})

test_that("HDRDA's statistics match manual values when (lambda, gamma) = (0.5, 0.75) and convex", {
  set.seed(42)
  p <- 250
  n1 <- 10
  n2 <- 15
  n3 <- 25

  iris_x <- as.matrix(iris[, -5])

  # The first cov matrix is the identity matrix.
  # To generate sigma2 and sigma3, we grab eigenvectors from the iris data set.
  # The eigenvalues are randomly generated.
  eigen2 <- eigen(tcrossprod(iris_x[rep(1:100, length=p), ]), symmetric=TRUE)
  eigen3 <- eigen(tcrossprod(iris_x[rep(30:130, length=p), ]), symmetric=TRUE)

  eigen2$values <- sort(rexp(p, rate=0.1), decreasing=TRUE)
  eigen3$values <- sort(rexp(p, rate=0.05), decreasing=TRUE)

  sigma2 <- with(eigen2, vectors %*% tcrossprod(diag(values), vectors))
  sigma3 <- with(eigen3, vectors %*% tcrossprod(diag(values), vectors))

  x1 <- rmvnorm(n1, mean=rep(0, p), sigma=diag(p))
  x2 <- rmvnorm(n2, mean=rep(3, p), sigma=sigma2)
  x3 <- rmvnorm(n3, mean=rep(-3, p), sigma=sigma3)

  x <- rbind(x1, x2, x3)
  y <- c(rep(1, n1), rep(2, n2), rep(3, n3))

  # Compute estimates for HDRDA manually.
  xbar1 <- colMeans(x1)
  xbar2 <- colMeans(x2)
  xbar3 <- colMeans(x3)

  x_centered1 <- scale(x1, center=TRUE, scale=FALSE)
  x_centered2 <- scale(x2, center=TRUE, scale=FALSE)
  x_centered3 <- scale(x3, center=TRUE, scale=FALSE)

  Sigma_eigen <- cov_eigen(x=x, y=y, pool=TRUE, fast=TRUE, tol=1e-6)
  D_q <- Sigma_eigen$values
  U1 <- Sigma_eigen$vectors
  q <- length(D_q)

  XU_1 <- x_centered1 %*% U1
  XU_2 <- x_centered2 %*% U1
  XU_3 <- x_centered3 %*% U1

  xbar1_U1 <- crossprod(U1, xbar1)
  xbar2_U1 <- crossprod(U1, xbar2)
  xbar3_U1 <- crossprod(U1, xbar3)

  Gamma <- 0.125 * diag(D_q) + 0.75 * diag(q)
  Gamma_inv <- diag(diag(Gamma)^(-1))

  W_1 <- 0.125 * cov_mle(XU_1) + Gamma
  W_2 <- 0.125 * cov_mle(XU_2) + Gamma
  W_3 <- 0.125 * cov_mle(XU_3) + Gamma

  W1_inv <- solve_chol(W_1)
  W2_inv <- solve_chol(W_2)
  W3_inv <- solve_chol(W_3)

  Q1 <- diag(n1) + 0.125 / n1 * XU_1 %*% tcrossprod(Gamma_inv, XU_1)
  Q2 <- diag(n2) + 0.125 / n2 * XU_2 %*% tcrossprod(Gamma_inv, XU_2)
  Q3 <- diag(n3) + 0.125 / n3 * XU_3 %*% tcrossprod(Gamma_inv, XU_3)

  hdrda_out <- rda_high_dim(x=x, y=y, lambda=0.5, gamma=0.75, shrinkage_type="convex")

  expect_equal(hdrda_out$q, q)
  expect_equal(hdrda_out$D_q, D_q)
  expect_equal(hdrda_out$U1, U1)

  expect_equal(hdrda_out$est[[1]]$n, n1)
  expect_equal(hdrda_out$est[[1]]$xbar, xbar1)
  expect_equal(hdrda_out$est[[1]]$prior, n1 / length(y))
  expect_equal(hdrda_out$est[[1]]$alpha, 0.25)
  expect_equal(hdrda_out$est[[1]]$XU, XU_1)
  expect_equal(hdrda_out$est[[1]]$xbar_U1, xbar1_U1)
  expect_equal(hdrda_out$est[[1]]$Gamma, diag(Gamma))
  expect_equal(hdrda_out$est[[1]]$Q, Q1)
  expect_equal(hdrda_out$est[[1]]$W_inv, W1_inv)

  expect_equal(hdrda_out$est[[2]]$n, n2)
  expect_equal(hdrda_out$est[[2]]$xbar, xbar2)
  expect_equal(hdrda_out$est[[2]]$prior, n2 / length(y))
  expect_equal(hdrda_out$est[[2]]$alpha, 0.25)
  expect_equal(hdrda_out$est[[2]]$XU, XU_2)
  expect_equal(hdrda_out$est[[2]]$xbar_U1, xbar2_U1)
  expect_equal(hdrda_out$est[[2]]$Gamma, diag(Gamma))
  expect_equal(hdrda_out$est[[2]]$Q, Q2)
  expect_equal(hdrda_out$est[[2]]$W_inv, W2_inv)

  expect_equal(hdrda_out$est[[3]]$n, n3)
  expect_equal(hdrda_out$est[[3]]$xbar, xbar3)
  expect_equal(hdrda_out$est[[3]]$prior, n3 / length(y))
  expect_equal(hdrda_out$est[[3]]$alpha, 0.25)
  expect_equal(hdrda_out$est[[3]]$XU, XU_3)
  expect_equal(hdrda_out$est[[3]]$xbar_U1, xbar3_U1)
  expect_equal(hdrda_out$est[[3]]$Gamma, diag(Gamma))
  expect_equal(hdrda_out$est[[3]]$Q, Q3)
  expect_equal(hdrda_out$est[[3]]$W_inv, W3_inv)
})

test_that("HDRDA posterior probabilities sum to one. (Issue #34)", {
  trn <- two_class_sim_data[1:100,]
  tst <- two_class_sim_data[101:105,]

  mod <- rda_high_dim(x = as.matrix(trn[, -ncol(trn)]), y = trn$Class)

  predict_out <- predict(mod, newdata=tst[, -ncol(tst)])
  posterior_probs <- predict_out$posterior
  scores <- predict_out$scores

  ones <- rep(1, nrow(tst))
  expect_equal(as.vector(rowSums(posterior_probs)), ones)
  expect_equal(colnames(posterior_probs), levels(two_class_sim_data$Class))
  expect_equal(colnames(scores), levels(two_class_sim_data$Class))
})

test_that("HDRDA correctly predicts one observation. (Issue #34)", {
  trn <- two_class_sim_data[1:100,]
  tst <- two_class_sim_data[101,]

  mod <- rda_high_dim(x = as.matrix(trn[, -ncol(trn)]), y = trn$Class)

  predict_out <- predict(mod, newdata=tst[, -ncol(tst)])
  posterior_probs <- predict_out$posterior
  scores <- predict_out$scores

  expect_equal(sum(posterior_probs), 1)
  expect_equal(names(posterior_probs), levels(two_class_sim_data$Class))
  expect_equal(names(scores), levels(two_class_sim_data$Class))
})

# Related to issue #41
test_that("The HDRDA classifier works properly when 1 feature used", {
  require('MASS')

  set.seed(42)
  n <- nrow(iris)
  train <- sample(seq_len(n), n / 2)
  n_test <- n - length(train)

  hdrda_out <- rda_high_dim(x = iris[train, 1], y = iris[train, 5], lambda=0.5, gamma=0.5)
  predicted <- predict(hdrda_out, iris[-train, 1])

  expect_equal(length(predicted$class), n_test)
  expect_is(predicted$posterior, "matrix")
  expect_equal(dim(predicted$posterior), c(n_test, nlevels(iris$Species)))
})
