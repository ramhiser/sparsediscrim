library('diagdiscrim')
library('testthat')

context("Simultaneous Diagonalization with Generalized Eigenvalues")

# TODO: Add this unit test back in when we generalize to K >= 2.
#test_that("more than 2 classes in data raises errors", {
#	num.classes <- 3
#	N <- 5
#	p <- 3
#	X <- rmvnorm(num.classes * N, mean = rep(0, p))
#	df <- data.frame(labels = sort(rep(seq_len(num.classes), n)), X)
#	expect_that(geneigen.cov(df), throws_error("K == 2 is not TRUE"))
#})

test_that("Two Low-Dimensional Covariance Matrices are Simultaneously Diagonalized", {
	set.seed(42)
	p <- 10
	x1 <- replicate(p, rnorm(20))
	x2 <- replicate(p, rnorm(20, 2))
	cov_1 <- cov(x1)
	cov_2 <- cov(x2)
	
	cov_2_eigen <- eigen(cov_2, symmetric = TRUE)
	Q_B <- cov_2_eigen$vectors %*% diag(cov_2_eigen$values^(-1/2))

	# We introduce an additional step to store the eigenvalues as well
	# for unit testing.
	#Q <- Q_B %*% eigen(t(Q_B) %*% cov_1 %*% Q_B, symmetric = TRUE)$vectors
	Q_A_eigen <- eigen(t(Q_B) %*% cov_1 %*% Q_B, symmetric = TRUE)
	Q <- Q_B %*% Q_A_eigen$vectors
	expect_equal(t(Q) %*% cov_1 %*% Q, diag(Q_A_eigen$values))
	expect_equal(t(Q) %*% cov_2 %*% Q, diag(p))
	
	# Cov(X_1 Q) should equal both the identity and Q' \Sigma_1 Q.
	expect_equal(
		cov(x1 %*% Q),
		diag(Q_A_eigen$values)
	)
	# Cov(X_2 Q) should equal both the identity and Q' \Sigma_2 Q = I_p.
	expect_equal(
		cov(x2 %*% Q),
		diag(p)
	)
})

test_that("Two High-Dimensional Covariance Matrices are Simultaneously Diagonalized", {
	set.seed(42)
	p <- 100
	x1 <- replicate(p, rnorm(20))
	x2 <- replicate(p, rnorm(20, 2))
	cov_1 <- cov(x1)
	cov_2 <- cov(x2)
	
	cov_2_eigen <- eigen(cov_2, symmetric = TRUE)
	
	# We note that cov_2_eigen is singular here and denote its rank
	# as q, where 1 <= q < p.
	# We apply our dimension reduction method that uses the top
	# k eigenvalues, where 1 <= k <= q. One approach to determining
	# k is to set it equal to the number of nonzero eigenvalues.
	# We apply the square root to R's built-in non-zero threshold
	# to declare an eigenvalue as 0.
	k <- sum(cov_2_eigen$values > sqrt(.Machine$double.eps))
	Q_B <- cov_2_eigen$vectors[,1:k] %*% diag(cov_2_eigen$values[1:k]^(-1/2))
	
	# The matrix Q_B should be have dimensions p x k.
	expect_equal(nrow(Q_B), p)
	expect_equal(ncol(Q_B), k)
	expect_equal(t(Q_B) %*% cov_2 %*% Q_B, diag(k))

	# We introduce an additional step to store the eigenvalues as well
	# for unit testing.
	#Q <- Q_B %*% eigen(t(Q_B) %*% cov_1 %*% Q_B, symmetric = TRUE)$vectors
	Q_A_eigen <- eigen(t(Q_B) %*% cov_1 %*% Q_B, symmetric = TRUE)
	Q <- Q_B %*% Q_A_eigen$vectors
	expect_equal(t(Q) %*% cov_1 %*% Q, diag(Q_A_eigen$values))
	expect_equal(t(Q) %*% cov_2 %*% Q, diag(k))
	
	# Cov(X_1 Q) should equal both the identity and Q' \Sigma_1 Q.
	expect_equal(
		cov(x1 %*% Q),
		diag(Q_A_eigen$values)
	)
	# Cov(X_2 Q) should equal both the identity and Q' \Sigma_2 Q = I_p.
	expect_equal(
		cov(x2 %*% Q),
		diag(k)
	)
})


test_that("The Function simdiag Simultaneously Diagonalizes Two Low-Dimensional Covariance Matrices", {
	set.seed(42)
	p <- 10
	x1 <- replicate(p, rnorm(20))
	x2 <- replicate(p, rnorm(20, 2))
	cov_1 <- cov(x1)
	cov_2 <- cov(x2)
	
	Q <- simdiag(cov_1, cov_2)$Q
	
	cov_2_eigen <- eigen(cov_2, symmetric = TRUE)
	Q_B <- cov_2_eigen$vectors %*% diag(cov_2_eigen$values^(-1/2))
	Q_A_eigen <- eigen(t(Q_B) %*% cov_1 %*% Q_B, symmetric = TRUE)

	expect_equal(t(Q) %*% cov_1 %*% Q, diag(Q_A_eigen$values))
	expect_equal(t(Q) %*% cov_2 %*% Q, diag(p))
	
	# Cov(X_1 Q) should equal both the identity and Q' \Sigma_1 Q.
	expect_equal(
		cov(x1 %*% Q),
		diag(Q_A_eigen$values)
	)
	# Cov(X_2 Q) should equal both the identity and Q' \Sigma_2 Q = I_p.
	expect_equal(
		cov(x2 %*% Q),
		diag(p)
	)
})

test_that("The Function simdiag Simultaneously Diagonalizes Two High-Dimensional Covariance Matrices", {
	set.seed(42)
	p <- 100
	x1 <- replicate(p, rnorm(20))
	x2 <- replicate(p, rnorm(20, 2))
	cov_1 <- cov(x1)
	cov_2 <- cov(x2)
	
	# Without a reduced dimension, the simdiag function should throw an error.
	expect_error(Q <- simdiag(cov_1, cov_2)$Q)
	expect_error(Q <- simdiag(cov_1, cov_2, dim_reduce = F)$Q)
	
	simdiag_out <- simdiag(cov_1, cov_2, dim_reduce = T)
	
	cov_2_eigen <- eigen(cov_2, symmetric = TRUE)
	k <- sum(cov_2_eigen$values > sqrt(.Machine$double.eps))
	expect_equal(simdiag_out$k, k)
	
	Q_B <- cov_2_eigen$vectors[,seq_len(k)] %*% diag(cov_2_eigen$values[seq_len(k)]^(-1/2))
	Q_A_eigen <- eigen(t(Q_B) %*% cov_1 %*% Q_B, symmetric = TRUE)
	Q <- Q_B %*% Q_A_eigen$vectors

	expect_equal(simdiag_out$Q, Q)
	expect_equal(t(simdiag_out$Q) %*% cov_1 %*% simdiag_out$Q, diag(Q_A_eigen$values[seq_len(k)]))
	expect_equal(t(simdiag_out$Q) %*% cov_2 %*% simdiag_out$Q, diag(k))
	
	# Cov(X_1 Q) should equal both the identity and Q' \Sigma_1 Q.
	expect_equal(
		cov(x1 %*% Q),
		diag(Q_A_eigen$values)
	)
	# Cov(X_2 Q) should equal both the identity and Q' \Sigma_2 Q = I_p.
	expect_equal(
		cov(x2 %*% Q),
		diag(k)
	)

	# Now, we repeat the above by manually selecting k to be 2.
	k <- 2
	simdiag_out <- simdiag(cov_1, cov_2, dim_reduce = T, k = k)
	cov_2_eigen <- eigen(cov_2, symmetric = TRUE)
	expect_equal(simdiag_out$k, k)
	
	Q_B <- cov_2_eigen$vectors[,seq_len(k)] %*% diag(cov_2_eigen$values[seq_len(k)]^(-1/2))
	Q_A_eigen <- eigen(t(Q_B) %*% cov_1 %*% Q_B, symmetric = TRUE)
	Q <- Q_B %*% Q_A_eigen$vectors

	expect_equal(simdiag_out$Q, Q)
	expect_equal(t(simdiag_out$Q) %*% cov_1 %*% simdiag_out$Q, diag(Q_A_eigen$values[seq_len(k)]))
	expect_equal(t(simdiag_out$Q) %*% cov_2 %*% simdiag_out$Q, diag(k))
	
	# Cov(X_1 Q) should equal both the identity and Q' \Sigma_1 Q.
	expect_equal(
		cov(x1 %*% Q),
		diag(Q_A_eigen$values)
	)
	# Cov(X_2 Q) should equal both the identity and Q' \Sigma_2 Q = I_p.
	expect_equal(
		cov(x2 %*% Q),
		diag(k)
	)
})

# There were several errors when k = 1 (either directly as specified by user or
#	indirectly when the eigen_pct or rank yielded 1.)
test_that("The function simdiag works correctly with k = 1", {
	set.seed(42)
	p <- 100
	x1 <- replicate(p, rnorm(20))
	x2 <- replicate(p, rnorm(20, 2))
	cov_1 <- cov(x1)
	cov_2 <- cov(x2)
		
	# Now, we repeat the above by manually selecting k to be 1.
	k <- 1
	simdiag_out <- simdiag(cov_1, cov_2, dim_reduce = T, k = k)
	cov_2_eigen <- eigen(cov_2, symmetric = TRUE)
	expect_equal(simdiag_out$k, k)
	
	Q_B <- cov_2_eigen$vectors[,1] * cov_2_eigen$values[1]^(-1/2)
	Q_A_eigen <- eigen(t(Q_B) %*% cov_1 %*% Q_B, symmetric = TRUE)
	Q <- Q_B %*% Q_A_eigen$vectors

	expect_equal(simdiag_out$Q, Q)
	expect_equal(drop(t(simdiag_out$Q) %*% cov_1 %*% simdiag_out$Q), Q_A_eigen$values[seq_len(k)])
	expect_equal(t(simdiag_out$Q) %*% cov_2 %*% simdiag_out$Q, diag(k))
	
	# Cov(X_1 Q) should equal both the identity and Q' \Sigma_1 Q.
	expect_equal(
		drop(cov(x1 %*% Q)),
		Q_A_eigen$values
	)
	# Cov(X_2 Q) should equal both the identity and Q' \Sigma_2 Q = I_p.
	expect_equal(
		cov(x2 %*% Q),
		diag(k)
	)
})