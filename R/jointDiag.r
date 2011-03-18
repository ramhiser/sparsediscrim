jointDiag.cov <- function(df, shrink = FALSE, shrink.val = 0.01) {
	if(shrink == FALSE) {
		shrink.val <- 0
	}
	covs <- daply(df, .(labels), function(df_k) {
		n_k <- nrow(df_k)
		p <- ncol(df_k) - 1
		(1 - shrink.val) * (n_k - 1) * cov(df_k[,-1]) / n_k + shrink.val * diag(p)
	})
	covs <- aperm(covs, perm = c(2,3,1))

	dimnames(covs) <- NULL
	covs
}