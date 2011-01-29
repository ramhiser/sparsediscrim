eigendiag.cov <- function(df, shrink = FALSE, shrink.val = 0.01) {
	if(shrink == FALSE) {
		shrink.val <- 0
	}
	covs <- dlply(df, .(labels), function(class.df) {
		n.k <- nrow(class.df)
		p <- ncol(class.df) - 1
		cov.k <- (1 - shrink.val) * (n.k - 1) * cov(class.df[,-1]) / n.k + shrink.val * diag(p)
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