geneigen.cov <- function(df, shrink = FALSE, shrink.val = 0.01) {
	K <- nlevels(factor(df$labels))
	stopifnot(K == 2)	
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

# TODO: Implement geneigen via LAPACK. See solve and message from Peter Dalgaard
# See: http://tolstoy.newcastle.edu.au/R/help/05/06/6995.html
geneigen <- function(A, B) {
	eigen(solve(A) %*% B)
}

diagonalize.geneigen <- function(A, B) {
	generalized.eigen(A, B)$vectors
}