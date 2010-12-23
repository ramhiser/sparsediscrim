generalized.eigen.cov <- function(df, shrink = FALSE, shrink.val = 0.01) {
	K <- nlevels(factor(df$labels))
	stopifnot(K == 2)	
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

generalized.eigen <- function(A, B) {
	eigen(solve(A) %*% B)
}

diagonalize.general.eigen <- function(A, B) {
	generalized.eigen(A, B)$vectors
}