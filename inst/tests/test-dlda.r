library('diagdiscrim')
library('mvtnorm')
library('plyr')

context("Diagonalized Linear Discriminant Analysis")

test_that("DLDA computes correct estimators and predictions and predictions with 2 classes having identity covariance matrix and p = 3", {
	set.seed(42)
	n1 <- n2 <- 10
	p <- 3
	mu1 <- 0
	mu2 <- 2

	x1 <- rmvnorm(n1, rep(mu1, p))
	x2 <- rmvnorm(n2, rep(mu2, p))
	
  x <- rbind(x1, x2)
  y <- gl(2, n1)
	
	xbar1 <- as.vector(aaply(x1, 2, mean))
	xbar2 <- as.vector(aaply(x2, 2, mean))
	
	cov1 <- (n1 - 1) * cov(x1) / n1
	cov2 <- (n2 - 1) * cov(x2) / n2
	
	cov_pool <- (n1 * cov1 + n2 * cov2) / (n1 + n2)
	
	dlda_out <- dlda(x = x, y = y)
	
	expect_that(all.equal(as.vector(dlda_out$var_pool), diag(cov_pool)), is_true())
	expect_that(all.equal(as.vector(dlda_out$var_pool), diag(cov_pool)), is_true())
	
	expect_that(all.equal(as.vector(dlda_out$est$`1`$xbar), xbar1), is_true())
	expect_that(all.equal(as.vector(dlda_out$est$`2`$xbar), xbar2), is_true())
	
	test1 <- test2 <- 500
	test_y <- c(rep(1, test1), rep(2, test2))
	test_x1 <- rmvnorm(test1, rep(mu1, p))
	test_x2 <- rmvnorm(test2, rep(mu2, p))
	test_x <- rbind(test_x1, test_x2)
	predictions <- predict(dlda_out, test_x)$class
	test_error <- mean(predictions != test_y)
	
	expect_that(test_error, equals(0.03))
})

test_that("DLDA computes correct estimators and predictions with 2 classes having identity covariance matrix and p = 20", {
	set.seed(42)
	n1 <- n2 <- 15
	p <- 20
	mu1 <- 0
	mu2 <- 2

	x1 <- rmvnorm(n1, rep(mu1, p))
	x2 <- rmvnorm(n2, rep(mu2, p))
	
  x <- rbind(x1, x2)
  y <- gl(2, n1)
	
	xbar1 <- as.vector(aaply(x1, 2, mean))
	xbar2 <- as.vector(aaply(x2, 2, mean))
	
	cov1 <- (n1 - 1) * cov(x1) / n1
	cov2 <- (n2 - 1) * cov(x2) / n2
	
	cov_pool <- (n1 * cov1 + n2 * cov2) / (n1 + n2)
	
	dlda_out <- dlda(x = x, y = y)
	
	expect_that(all.equal(as.vector(dlda_out$var_pool), diag(cov_pool)), is_true())
	expect_that(all.equal(as.vector(dlda_out$var_pool), diag(cov_pool)), is_true())
	
	expect_that(all.equal(as.vector(dlda_out$est$`1`$xbar), xbar1), is_true())
	expect_that(all.equal(as.vector(dlda_out$est$`2`$xbar), xbar2), is_true())
	
	test1 <- test2 <- 500
	test_y <- c(rep(1, test1), rep(2, test2))
	test_x1 <- rmvnorm(test1, rep(mu1, p))
	test_x2 <- rmvnorm(test2, rep(mu2, p))
	test_x <- rbind(test_x1, test_x2)
	predictions <- predict(dlda_out, test_x)$class
	test_error <- mean(predictions != test_y)
	
	expect_that(test_error, equals(0))
})

test_that("DLDA computes correct estimators and predictions with 3 classes having identity covariance matrix and p = 5", {
	set.seed(42)
	n1 <- n2 <- n3 <- 15
	p <- 5
	mu1 <- 0
	mu2 <- 2
	mu3 <- 4

	x1 <- rmvnorm(n1, rep(mu1, p))
	x2 <- rmvnorm(n2, rep(mu2, p))
	x3 <- rmvnorm(n3, rep(mu3, p))
	
  x <- rbind(x1, x2, x3)
  y <- gl(3, n1)
	
	xbar1 <- as.vector(aaply(x1, 2, mean))
	xbar2 <- as.vector(aaply(x2, 2, mean))
	xbar3 <- as.vector(aaply(x3, 2, mean))
	
	cov1 <- (n1 - 1) * cov(x1) / n1
	cov2 <- (n2 - 1) * cov(x2) / n2
	cov3 <- (n3 - 1) * cov(x3) / n3
	
	cov_pool <- (n1 * cov1 + n2 * cov2 + n3 * cov3) / (n1 + n2 + n3)
	
	dlda_out <- dlda(x = x, y = y)
	
	expect_that(all.equal(as.vector(dlda_out$var_pool), diag(cov_pool)), is_true())
	expect_that(all.equal(as.vector(dlda_out$var_pool), diag(cov_pool)), is_true())
	expect_that(all.equal(as.vector(dlda_out$var_pool), diag(cov_pool)), is_true())
	
	expect_that(all.equal(as.vector(dlda_out$est$`1`$xbar), xbar1), is_true())
	expect_that(all.equal(as.vector(dlda_out$est$`2`$xbar), xbar2), is_true())
	expect_that(all.equal(as.vector(dlda_out$est$`3`$xbar), xbar3), is_true())
	
	test1 <- test2 <- test3 <- 500
	test_y <- c(rep(1, test1), rep(2, test2), rep(3, test3))
	test_x1 <- rmvnorm(test1, rep(mu1, p))
	test_x2 <- rmvnorm(test2, rep(mu2, p))
	test_x3 <- rmvnorm(test3, rep(mu3, p))
	test_x <- rbind(test_x1, test_x2, test_x3)
	predictions <- predict(dlda_out, test_x)$class
	test_error <- mean(predictions != test_y)
	
	expect_that(test_error, equals(0.0226666667))
})

test_that("DLDA computes correct estimators and predictions with 3 classes having identity covariance matrix and p = 50", {
	set.seed(42)
	n1 <- n2 <- n3 <- 15
	p <- 50
	mu1 <- 0
	mu2 <- 2
	mu3 <- 4

	x1 <- rmvnorm(n1, rep(mu1, p))
	x2 <- rmvnorm(n2, rep(mu2, p))
	x3 <- rmvnorm(n3, rep(mu3, p))
	
  x <- rbind(x1, x2, x3)
  y <- gl(3, n1)

	xbar1 <- as.vector(aaply(x1, 2, mean))
	xbar2 <- as.vector(aaply(x2, 2, mean))
	xbar3 <- as.vector(aaply(x3, 2, mean))
	
	cov1 <- (n1 - 1) * cov(x1) / n1
	cov2 <- (n2 - 1) * cov(x2) / n2
	cov3 <- (n3 - 1) * cov(x3) / n3
	
	cov_pool <- (n1 * cov1 + n2 * cov2 + n3 * cov3) / (n1 + n2 + n3)
	
	dlda_out <- dlda(x = x, y = y)
	
	expect_that(all.equal(as.vector(dlda_out$var_pool), diag(cov_pool)), is_true())
	expect_that(all.equal(as.vector(dlda_out$var_pool), diag(cov_pool)), is_true())
	expect_that(all.equal(as.vector(dlda_out$var_pool), diag(cov_pool)), is_true())
	
	expect_that(all.equal(as.vector(dlda_out$est$`1`$xbar), xbar1), is_true())
	expect_that(all.equal(as.vector(dlda_out$est$`2`$xbar), xbar2), is_true())
	expect_that(all.equal(as.vector(dlda_out$est$`3`$xbar), xbar3), is_true())
	
	test1 <- test2 <- test3 <- 500
	test_y <- c(rep(1, test1), rep(2, test2), rep(3, test3))
	test_x1 <- rmvnorm(test1, rep(mu1, p))
	test_x2 <- rmvnorm(test2, rep(mu2, p))
	test_x3 <- rmvnorm(test3, rep(mu3, p))
	test_x <- rbind(test_x1, test_x2, test_x3)
	predictions <- predict(dlda_out, test_x)$class
	test_error <- mean(predictions != test_y)
	
	expect_that(test_error, equals(0))
})

test_that("DLDA computes correct estimators and predictions with 5 classes having intraclass covariance matrix and p = 5", {
	set.seed(42)
	n1 <- n2 <- n3 <- n4 <- n5 <- 50
	p <- 5
	mu1 <- 0
	mu2 <- 2
	mu3 <- 4
	mu4 <- 6
	mu5 <- 8
	
	rho <- 0.9
	
	Sigma <- rho * matrix(1, nrow = p, ncol = p) + (1 - rho) * diag(p)

	x1 <- rmvnorm(n1, rep(mu1, p), sigma = Sigma)
	x2 <- rmvnorm(n2, rep(mu2, p), sigma = Sigma)
	x3 <- rmvnorm(n3, rep(mu3, p), sigma = Sigma)
	x4 <- rmvnorm(n4, rep(mu4, p), sigma = Sigma)
	x5 <- rmvnorm(n5, rep(mu5, p), sigma = Sigma)
	
  x <- rbind(x1, x2, x3, x4, x5)
  y <- gl(5, n1)

	xbar1 <- as.vector(aaply(x1, 2, mean))
	xbar2 <- as.vector(aaply(x2, 2, mean))
	xbar3 <- as.vector(aaply(x3, 2, mean))
	xbar4 <- as.vector(aaply(x4, 2, mean))
	xbar5 <- as.vector(aaply(x5, 2, mean))
	
	cov1 <- (n1 - 1) * cov(x1) / n1
	cov2 <- (n2 - 1) * cov(x2) / n2
	cov3 <- (n3 - 1) * cov(x3) / n3
	cov4 <- (n4 - 1) * cov(x4) / n4
	cov5 <- (n5 - 1) * cov(x5) / n5
	
	cov_pool <- (n1 * cov1 + n2 * cov2 + n3 * cov3 + n4 * cov4 + n5 * cov5) / (n1 + n2 + n3 + n4 + n5)
	
	dlda_out <- dlda(x = x, y = y)
	
	expect_that(all.equal(as.vector(dlda_out$var_pool), diag(cov_pool)), is_true())
	expect_that(all.equal(as.vector(dlda_out$var_pool), diag(cov_pool)), is_true())
	expect_that(all.equal(as.vector(dlda_out$var_pool), diag(cov_pool)), is_true())
	expect_that(all.equal(as.vector(dlda_out$var_pool), diag(cov_pool)), is_true())
	expect_that(all.equal(as.vector(dlda_out$var_pool), diag(cov_pool)), is_true())
	
	expect_that(all.equal(as.vector(dlda_out$est$`1`$xbar), xbar1), is_true())
	expect_that(all.equal(as.vector(dlda_out$est$`2`$xbar), xbar2), is_true())
	expect_that(all.equal(as.vector(dlda_out$est$`3`$xbar), xbar3), is_true())
	expect_that(all.equal(as.vector(dlda_out$est$`4`$xbar), xbar4), is_true())
	expect_that(all.equal(as.vector(dlda_out$est$`5`$xbar), xbar5), is_true())
	
	test1 <- test2 <- test3 <- test4 <- test5 <- 500
	test_y <- c(rep(1, test1), rep(2, test2), rep(3, test3), rep(4, test4), rep(5, test5))
	test_x1 <- rmvnorm(test1, rep(mu1, p))
	test_x2 <- rmvnorm(test2, rep(mu2, p))
	test_x3 <- rmvnorm(test3, rep(mu3, p))
	test_x4 <- rmvnorm(test4, rep(mu4, p))
	test_x5 <- rmvnorm(test5, rep(mu5, p))
	test_x <- rbind(test_x1, test_x2, test_x3, test_x4, test_x5)
	predictions <- predict(dlda_out, test_x)$class
	test_error <- mean(predictions != test_y)
	
	expect_that(test_error, equals(0.022))
})

test_that("DLDA computes correct estimators and predictions with 5 classes having intraclass covariance matrix and p = 50", {
	set.seed(42)
	n1 <- n2 <- n3 <- n4 <- n5 <- 50
	p <- 50
	mu1 <- 0
	mu2 <- 2
	mu3 <- 4
	mu4 <- 6
	mu5 <- 8
	
	rho <- 0.9
	
	Sigma <- rho * matrix(1, nrow = p, ncol = p) + (1 - rho) * diag(p)

	x1 <- rmvnorm(n1, rep(mu1, p), sigma = Sigma)
	x2 <- rmvnorm(n2, rep(mu2, p), sigma = Sigma)
	x3 <- rmvnorm(n3, rep(mu3, p), sigma = Sigma)
	x4 <- rmvnorm(n4, rep(mu4, p), sigma = Sigma)
	x5 <- rmvnorm(n5, rep(mu5, p), sigma = Sigma)
	
  x <- rbind(x1, x2, x3, x4, x5)
  y <- gl(5, n1)

	xbar1 <- as.vector(aaply(x1, 2, mean))
	xbar2 <- as.vector(aaply(x2, 2, mean))
	xbar3 <- as.vector(aaply(x3, 2, mean))
	xbar4 <- as.vector(aaply(x4, 2, mean))
	xbar5 <- as.vector(aaply(x5, 2, mean))
	
	cov1 <- (n1 - 1) * cov(x1) / n1
	cov2 <- (n2 - 1) * cov(x2) / n2
	cov3 <- (n3 - 1) * cov(x3) / n3
	cov4 <- (n4 - 1) * cov(x4) / n4
	cov5 <- (n5 - 1) * cov(x5) / n5
	
	cov_pool <- (n1 * cov1 + n2 * cov2 + n3 * cov3 + n4 * cov4 + n5 * cov5) / (n1 + n2 + n3 + n4 + n5)
	
	dlda_out <- dlda(x = x, y = y)
	
	expect_that(all.equal(as.vector(dlda_out$var_pool), diag(cov_pool)), is_true())
	expect_that(all.equal(as.vector(dlda_out$var_pool), diag(cov_pool)), is_true())
	expect_that(all.equal(as.vector(dlda_out$var_pool), diag(cov_pool)), is_true())
	expect_that(all.equal(as.vector(dlda_out$var_pool), diag(cov_pool)), is_true())
	expect_that(all.equal(as.vector(dlda_out$var_pool), diag(cov_pool)), is_true())
	
	expect_that(all.equal(as.vector(dlda_out$est$`1`$xbar), xbar1), is_true())
	expect_that(all.equal(as.vector(dlda_out$est$`2`$xbar), xbar2), is_true())
	expect_that(all.equal(as.vector(dlda_out$est$`3`$xbar), xbar3), is_true())
	expect_that(all.equal(as.vector(dlda_out$est$`4`$xbar), xbar4), is_true())
	expect_that(all.equal(as.vector(dlda_out$est$`5`$xbar), xbar5), is_true())
	
	test1 <- test2 <- test3 <- test4 <- test5 <- 500
	test_y <- c(rep(1, test1), rep(2, test2), rep(3, test3), rep(4, test4), rep(5, test5))
	test_x1 <- rmvnorm(test1, rep(mu1, p))
	test_x2 <- rmvnorm(test2, rep(mu2, p))
	test_x3 <- rmvnorm(test3, rep(mu3, p))
	test_x4 <- rmvnorm(test4, rep(mu4, p))
	test_x5 <- rmvnorm(test5, rep(mu5, p))
	test_x <- rbind(test_x1, test_x2, test_x3, test_x4, test_x5)
	predictions <- predict(dlda_out, test_x)$class
	test_error <- mean(predictions != test_y)
	
	expect_that(test_error, equals(0))
})





