eigendiag.cov <- function(df, shrink = FALSE, shrink.val = 0.01) {
	if(shrink == FALSE) {
		shrink.val <- 0
	}
	covs <- dlply(df, .(labels), function(df_k) {
		n_k <- nrow(df_k)
		p <- ncol(df_k) - 1
		cov.k <- (1 - shrink.val) * (n_k - 1) * cov(df_k[,-1]) / n_k + shrink.val * diag(p)
		dimnames(cov.k) <- NULL
		cov.k
	})
	covs
}

eigendiag <- function(mat.list) {
	eigen.list <- lapply(sig.list, function(x) eigen(x)$vectors)
	Q <- Reduce("+", eigen.list)
	list(Q = Q, eigenvectors.list = eigen.list)
}