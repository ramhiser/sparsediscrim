pool.cov <- function(df, shrink = FALSE, shrink.val = 0.01) {
	if(shrink == FALSE) {
		shrink.val <- 0
	}
	N <- nrow(df)
	p <- ncol(df_k) - 1
	covs <- dlply(df, .(labels), function(df_k) {
		n_k <- nrow(df_k)
		cov.k <- (n_k - 1) * data.matrix(cov(df_k[,-1]))
		cov.k
	})
	cov_pool <- (1 - shrink.val) * Reduce("+", covs) / N + shrink.val * diag(p)
	covs
}

diagonalize.pool.cov <- function(cov_pooled, eigen_tol = 1e-6) {
	cov_pooled_eigen <- eigen(cov_pooled, symmetric = TRUE)
	t(cov_pooled_eigen$vectors)
}