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
	p <- 10
	x1 <- replicate(p, rnorm(20))
	x2 <- replicate(p, rnorm(20, 2))
	cov_1 <- cov(x1)
	cov_2 <- cov(x2)
	T_chol <- chol(cov_1)
	
	# (T')^{-1} \Sigma_1 T^{-1} should equal the identity matrix.
	expect_equal(
		solve(t(T_chol)) %*% cov_1 %*% solve(T_chol),
		diag(p)
	)
	
	E <- solve(T_chol)
	B_star <- t(E) %*% cov_2 %*% E
	C <- eigen(B_star, symmetric = T)$vectors
	F_mat <- E %*% C
	
	# F' \Sigma_1 F should equal the identity matrix.
	expect_equal(
		t(F_mat) %*% cov_1 %*% F_mat,
		diag(p)
	)

	# F' \Sigma_2 F should equal the diagonal matrix with generalized evals
	# of Bx = lambda Ax. Here, A = cov_1 and B = cov_2.
	# For the case where A is positive definite, we can find generalized
	# eigenvalues by finding the eigenvalues of A^{-1} B x = lambda x.
	expect_equal(
		t(F_mat) %*% cov_2 %*% F_mat,
		diag(eigen(solve(cov_1) %*% cov_2)$values)
	)
	
	# Cov(X_1 F) should equal both the identity and F' \Sigma_1 F.
	expect_equal(
		cov(x1 %*% F_mat),
		diag(p)
	)
	expect_equal(
		cov(x1 %*% F_mat),
		t(F_mat) %*% cov_1 %*% F_mat
	)
	
	# Cov(X_2 F) should equal both the matrix of eigenvalues of A^{-1} B and F' \Sigma_2 F.
	expect_equal(
		cov(x2 %*% F_mat),
		diag(eigen(solve(cov_1) %*% cov_2)$values)
	)
	expect_equal(
		cov(x2 %*% F_mat),
		t(F_mat) %*% cov_2 %*% F_mat
	)
})

test_that("Two High-Dimensional Covariance Matrices are Simultaneously Diagonalized", {
	p <- 100
	x1 <- replicate(p, rnorm(20))
	x2 <- replicate(p, rnorm(20, 2))
	cov_1 <- cov(x1)
	cov_2 <- cov(x2)
	cov_1 <- cov_1 + .01 * diag(p)
	T_chol <- chol(cov_1)
	
	# (T')^{-1} \Sigma_1 T^{-1} should equal the identity matrix.
	expect_equal(
		solve(t(T_chol)) %*% cov_1 %*% solve(T_chol),
		diag(p)
	)
	
	E <- solve(T_chol)
	B_star <- t(E) %*% cov_2 %*% E
	C <- eigen(B_star, symmetric = T)$vectors
	F_mat <- E %*% C
	
	# F' \Sigma_1 F should equal the identity matrix.
	expect_equal(
		t(F_mat) %*% cov_1 %*% F_mat,
		diag(p)
	)

	# F' \Sigma_2 F should equal the diagonal matrix with generalized evals
	# of Bx = lambda Ax. Here, A = cov_1 and B = cov_2.
	# For the case where A is positive definite, we can find generalized
	# eigenvalues by finding the eigenvalues of A^{-1} B x = lambda x.
	expect_equal(
		t(F_mat) %*% cov_2 %*% F_mat,
		diag(eigen(solve(cov_1) %*% cov_2)$values)
	)
	
	# Cov(X_1 F) should equal both the identity and F' \Sigma_1 F.
	expect_equal(
		cov(x1 %*% F_mat),
		diag(p)
	)
	expect_equal(
		cov(x1 %*% F_mat),
		t(F_mat) %*% cov_1 %*% F_mat
	)
	
	# Cov(X_2 F) should equal both the matrix of eigenvalues of A^{-1} B and F' \Sigma_2 F.
	expect_equal(
		cov(x2 %*% F_mat),
		diag(eigen(solve(cov_1) %*% cov_2)$values)
	)
	expect_equal(
		cov(x2 %*% F_mat),
		t(F_mat) %*% cov_2 %*% F_mat
	)
})