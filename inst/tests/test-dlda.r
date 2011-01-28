library('diagdiscrim')
library('mvtnorm')

context("Diagonalized Linear Discriminant Analysis")

test_that("DLDA computes correct estimators and predictions and predictions with 2 classes having identity covariance matrix and p = 3", {
	set.seed(42)
	n1 <- n2 <- 10
	p <- 3
	mu1 <- 0
	mu2 <- 2

	x1 <- rmvnorm(n1, rep(mu1, p))
	x2 <- rmvnorm(n2, rep(mu2, p))
	
	train.df <- rbind.data.frame(cbind(1, x1), cbind(2, x2))
	names(train.df) <- c("labels", paste("X", seq_len(ncol(train.df) - 1), sep = ""))
	train.df$labels <- factor(train.df$labels)
	
	xbar1 <- as.vector(aaply(x1, 2, mean))
	xbar2 <- as.vector(aaply(x2, 2, mean))
	
	cov1 <- (n1 - 1) * cov(x1) / n1
	cov2 <- (n2 - 1) * cov(x2) / n2
	
	cov.pool <- (n1 * cov1 + n2 * cov2) / (n1 + n2)
	
	dlda.out <- dlda(train.df)
	
	expect_that(all.equal(as.vector(dlda.out$estimators$`1`$var), diag(cov.pool)), is_true())
	expect_that(all.equal(as.vector(dlda.out$estimators$`2`$var), diag(cov.pool)), is_true())
	
	expect_that(all.equal(as.vector(dlda.out$estimators$`1`$xbar), xbar1), is_true())
	expect_that(all.equal(as.vector(dlda.out$estimators$`2`$xbar), xbar2), is_true())
	
	test1 <- test2 <- 500
	test.labels <- c(rep(1, test1), rep(2, test2))
	test.x1 <- rmvnorm(test1, rep(mu1, p))
	test.x2 <- rmvnorm(test2, rep(mu2, p))
	test.x <- rbind(test.x1, test.x2)
	predictions <- predict(dlda.out, test.x)
	test.error <- mean(predictions != test.labels)
	
	expect_that(test.error, equals(0.03))
})

test_that("DLDA computes correct estimators and predictions with 2 classes having identity covariance matrix and p = 20", {
	set.seed(42)
	n1 <- n2 <- 15
	p <- 20
	mu1 <- 0
	mu2 <- 2

	x1 <- rmvnorm(n1, rep(mu1, p))
	x2 <- rmvnorm(n2, rep(mu2, p))
	
	train.df <- rbind.data.frame(cbind(1, x1), cbind(2, x2))
	names(train.df) <- c("labels", paste("X", seq_len(ncol(train.df) - 1), sep = ""))
	train.df$labels <- factor(train.df$labels)
	
	xbar1 <- as.vector(aaply(x1, 2, mean))
	xbar2 <- as.vector(aaply(x2, 2, mean))
	
	cov1 <- (n1 - 1) * cov(x1) / n1
	cov2 <- (n2 - 1) * cov(x2) / n2
	
	cov.pool <- (n1 * cov1 + n2 * cov2) / (n1 + n2)
	
	dlda.out <- dlda(train.df)
	
	expect_that(all.equal(as.vector(dlda.out$estimators$`1`$var), diag(cov.pool)), is_true())
	expect_that(all.equal(as.vector(dlda.out$estimators$`2`$var), diag(cov.pool)), is_true())
	
	expect_that(all.equal(as.vector(dlda.out$estimators$`1`$xbar), xbar1), is_true())
	expect_that(all.equal(as.vector(dlda.out$estimators$`2`$xbar), xbar2), is_true())
	
	test1 <- test2 <- 500
	test.labels <- c(rep(1, test1), rep(2, test2))
	test.x1 <- rmvnorm(test1, rep(mu1, p))
	test.x2 <- rmvnorm(test2, rep(mu2, p))
	test.x <- rbind(test.x1, test.x2)
	predictions <- predict(dlda.out, test.x)
	test.error <- mean(predictions != test.labels)
	
	expect_that(test.error, equals(0))
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
	
	train.df <- rbind.data.frame(cbind(1, x1), cbind(2, x2), cbind(3, x3))
	names(train.df) <- c("labels", paste("X", seq_len(ncol(train.df) - 1), sep = ""))
	train.df$labels <- factor(train.df$labels)
	
	xbar1 <- as.vector(aaply(x1, 2, mean))
	xbar2 <- as.vector(aaply(x2, 2, mean))
	xbar3 <- as.vector(aaply(x3, 2, mean))
	
	cov1 <- (n1 - 1) * cov(x1) / n1
	cov2 <- (n2 - 1) * cov(x2) / n2
	cov3 <- (n3 - 1) * cov(x3) / n3
	
	cov.pool <- (n1 * cov1 + n2 * cov2 + n3 * cov3) / (n1 + n2 + n3)
	
	dlda.out <- dlda(train.df)
	
	expect_that(all.equal(as.vector(dlda.out$estimators$`1`$var), diag(cov.pool)), is_true())
	expect_that(all.equal(as.vector(dlda.out$estimators$`2`$var), diag(cov.pool)), is_true())
	expect_that(all.equal(as.vector(dlda.out$estimators$`3`$var), diag(cov.pool)), is_true())
	
	expect_that(all.equal(as.vector(dlda.out$estimators$`1`$xbar), xbar1), is_true())
	expect_that(all.equal(as.vector(dlda.out$estimators$`2`$xbar), xbar2), is_true())
	expect_that(all.equal(as.vector(dlda.out$estimators$`3`$xbar), xbar3), is_true())
	
	test1 <- test2 <- test3 <- 500
	test.labels <- c(rep(1, test1), rep(2, test2), rep(3, test3))
	test.x1 <- rmvnorm(test1, rep(mu1, p))
	test.x2 <- rmvnorm(test2, rep(mu2, p))
	test.x3 <- rmvnorm(test3, rep(mu3, p))
	test.x <- rbind(test.x1, test.x2, test.x3)
	predictions <- predict(dlda.out, test.x)
	test.error <- mean(predictions != test.labels)
	
	expect_that(test.error, equals(0.0226666667))
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
	
	train.df <- rbind.data.frame(cbind(1, x1), cbind(2, x2), cbind(3, x3))
	names(train.df) <- c("labels", paste("X", seq_len(ncol(train.df) - 1), sep = ""))
	train.df$labels <- factor(train.df$labels)
	
	xbar1 <- as.vector(aaply(x1, 2, mean))
	xbar2 <- as.vector(aaply(x2, 2, mean))
	xbar3 <- as.vector(aaply(x3, 2, mean))
	
	cov1 <- (n1 - 1) * cov(x1) / n1
	cov2 <- (n2 - 1) * cov(x2) / n2
	cov3 <- (n3 - 1) * cov(x3) / n3
	
	cov.pool <- (n1 * cov1 + n2 * cov2 + n3 * cov3) / (n1 + n2 + n3)
	
	dlda.out <- dlda(train.df)
	
	expect_that(all.equal(as.vector(dlda.out$estimators$`1`$var), diag(cov.pool)), is_true())
	expect_that(all.equal(as.vector(dlda.out$estimators$`2`$var), diag(cov.pool)), is_true())
	expect_that(all.equal(as.vector(dlda.out$estimators$`3`$var), diag(cov.pool)), is_true())
	
	expect_that(all.equal(as.vector(dlda.out$estimators$`1`$xbar), xbar1), is_true())
	expect_that(all.equal(as.vector(dlda.out$estimators$`2`$xbar), xbar2), is_true())
	expect_that(all.equal(as.vector(dlda.out$estimators$`3`$xbar), xbar3), is_true())
	
	test1 <- test2 <- test3 <- 500
	test.labels <- c(rep(1, test1), rep(2, test2), rep(3, test3))
	test.x1 <- rmvnorm(test1, rep(mu1, p))
	test.x2 <- rmvnorm(test2, rep(mu2, p))
	test.x3 <- rmvnorm(test3, rep(mu3, p))
	test.x <- rbind(test.x1, test.x2, test.x3)
	predictions <- predict(dlda.out, test.x)
	test.error <- mean(predictions != test.labels)
	
	expect_that(test.error, equals(0))
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
	
	train.df <- rbind.data.frame(cbind(1, x1), cbind(2, x2), cbind(3, x3), cbind(4, x4), cbind(5, x5))
	names(train.df) <- c("labels", paste("X", seq_len(ncol(train.df) - 1), sep = ""))
	train.df$labels <- factor(train.df$labels)
	
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
	
	cov.pool <- (n1 * cov1 + n2 * cov2 + n3 * cov3 + n4 * cov4 + n5 * cov5) / (n1 + n2 + n3 + n4 + n5)
	
	dlda.out <- dlda(train.df)
	
	expect_that(all.equal(as.vector(dlda.out$estimators$`1`$var), diag(cov.pool)), is_true())
	expect_that(all.equal(as.vector(dlda.out$estimators$`2`$var), diag(cov.pool)), is_true())
	expect_that(all.equal(as.vector(dlda.out$estimators$`3`$var), diag(cov.pool)), is_true())
	expect_that(all.equal(as.vector(dlda.out$estimators$`4`$var), diag(cov.pool)), is_true())
	expect_that(all.equal(as.vector(dlda.out$estimators$`5`$var), diag(cov.pool)), is_true())
	
	expect_that(all.equal(as.vector(dlda.out$estimators$`1`$xbar), xbar1), is_true())
	expect_that(all.equal(as.vector(dlda.out$estimators$`2`$xbar), xbar2), is_true())
	expect_that(all.equal(as.vector(dlda.out$estimators$`3`$xbar), xbar3), is_true())
	expect_that(all.equal(as.vector(dlda.out$estimators$`4`$xbar), xbar4), is_true())
	expect_that(all.equal(as.vector(dlda.out$estimators$`5`$xbar), xbar5), is_true())
	
	test1 <- test2 <- test3 <- test4 <- test5 <- 500
	test.labels <- c(rep(1, test1), rep(2, test2), rep(3, test3), rep(4, test4), rep(5, test5))
	test.x1 <- rmvnorm(test1, rep(mu1, p))
	test.x2 <- rmvnorm(test2, rep(mu2, p))
	test.x3 <- rmvnorm(test3, rep(mu3, p))
	test.x4 <- rmvnorm(test4, rep(mu4, p))
	test.x5 <- rmvnorm(test5, rep(mu5, p))
	test.x <- rbind(test.x1, test.x2, test.x3, test.x4, test.x5)
	predictions <- predict(dlda.out, test.x)
	test.error <- mean(predictions != test.labels)
	
	expect_that(test.error, equals(0.022))
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
	
	train.df <- rbind.data.frame(cbind(1, x1), cbind(2, x2), cbind(3, x3), cbind(4, x4), cbind(5, x5))
	names(train.df) <- c("labels", paste("X", seq_len(ncol(train.df) - 1), sep = ""))
	train.df$labels <- factor(train.df$labels)
	
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
	
	cov.pool <- (n1 * cov1 + n2 * cov2 + n3 * cov3 + n4 * cov4 + n5 * cov5) / (n1 + n2 + n3 + n4 + n5)
	
	dlda.out <- dlda(train.df)
	
	expect_that(all.equal(as.vector(dlda.out$estimators$`1`$var), diag(cov.pool)), is_true())
	expect_that(all.equal(as.vector(dlda.out$estimators$`2`$var), diag(cov.pool)), is_true())
	expect_that(all.equal(as.vector(dlda.out$estimators$`3`$var), diag(cov.pool)), is_true())
	expect_that(all.equal(as.vector(dlda.out$estimators$`4`$var), diag(cov.pool)), is_true())
	expect_that(all.equal(as.vector(dlda.out$estimators$`5`$var), diag(cov.pool)), is_true())
	
	expect_that(all.equal(as.vector(dlda.out$estimators$`1`$xbar), xbar1), is_true())
	expect_that(all.equal(as.vector(dlda.out$estimators$`2`$xbar), xbar2), is_true())
	expect_that(all.equal(as.vector(dlda.out$estimators$`3`$xbar), xbar3), is_true())
	expect_that(all.equal(as.vector(dlda.out$estimators$`4`$xbar), xbar4), is_true())
	expect_that(all.equal(as.vector(dlda.out$estimators$`5`$xbar), xbar5), is_true())
	
	test1 <- test2 <- test3 <- test4 <- test5 <- 500
	test.labels <- c(rep(1, test1), rep(2, test2), rep(3, test3), rep(4, test4), rep(5, test5))
	test.x1 <- rmvnorm(test1, rep(mu1, p))
	test.x2 <- rmvnorm(test2, rep(mu2, p))
	test.x3 <- rmvnorm(test3, rep(mu3, p))
	test.x4 <- rmvnorm(test4, rep(mu4, p))
	test.x5 <- rmvnorm(test5, rep(mu5, p))
	test.x <- rbind(test.x1, test.x2, test.x3, test.x4, test.x5)
	predictions <- predict(dlda.out, test.x)
	test.error <- mean(predictions != test.labels)
	
	expect_that(test.error, equals(0))
})





