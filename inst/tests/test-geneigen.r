context("Simultaneous Diagonalization with Generalized Eigenvalues")

test_that("more than 2 classes in data raises errors", {
	num.classes <- 3
	N <- 5
	p <- 3
	X <- rmvnorm(num.classes * N, mean = rep(0, p))
	df <- data.frame(labels = sort(rep(seq_len(num.classes), n)), X)
	expect_that(geneigen.cov(df), throws_error("K == 2 is not TRUE"))
})

test_that("generalized eigenvalues of 2 predefined matrices equal precomputed values", {
	X1 <- matrix(c(1,0.5, 0.5, 2), nrow = 2, ncol = 2, byrow = TRUE)
	X2 <- matrix(c(2, -1, -1, 4), nrow = 2, ncol = 2, byrow = TRUE)

	B <- diagonalize.geneigen(X1, X2)

	expect_that(round(t(B) %*% X1 %*% B, digits = 3), equals(diag(c(0.862, 1.805))))
	expect_that(round(t(B) %*% X2 %*% B, digits = 3), equals(diag(c(3.609, 1.724))))
	
	set.seed(42)
	N <- 10
	p <- 5
	X1 <- crossprod(rmvnorm(N, mean = rep(0, p)))
	X2 <- crossprod(rmvnorm(N, mean = rep(1, p)))

	B <- diagonalize.geneigen(X1, X2)

	expect_that(round(t(B) %*% X1 %*% B, digits = 3), equals(diag(c(2.062, 11.883, 8.503, 5.442, 9.651))))
	expect_that(round(t(B) %*% X2 %*% B, digits = 3), equals(diag(c(47.788, 11.951, 5.924, 2.078, 0.964))))
})