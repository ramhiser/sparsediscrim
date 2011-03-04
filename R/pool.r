pool.cov <- function(df, shrink = FALSE, shrink.val = 0.01) {
	if(shrink == FALSE) {
		shrink.val <- 0
	}
	N <- nrow(df)
	p <- ncol(class.df) - 1
	covs <- dlply(df, .(labels), function(class.df) {
		n.k <- nrow(class.df)
		cov.k <- (n.k - 1) * data.matrix(cov(class.df[,-1]))
		cov.k
	})
	cov_pool <- (1 - shrink.val) * Reduce("+", covs) / N + shrink.val * diag(p)
	covs
}

diagonalize.pool.cov <- function(cov_pooled, eigen_tol = 1e-6) {
	cov_pooled_eigen <- eigen(cov_pooled, symmetric = TRUE)
	t(cov_pooled_eigen$vectors)
}