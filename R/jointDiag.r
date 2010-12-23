jointDiag.cov <- function(df, shrink = FALSE, shrink.val = 0.01) {
	if(shrink == FALSE) {
		shrink.val <- 0
	}
	covs <- daply(df, .(labels), function(class.df) {
		n.k <- nrow(class.df)
		p <- ncol(class.df) - 1
		(1 - shrink.val) * (n.k - 1) * cov(class.df[,-1]) / n.k + shrink.val * diag(p)
	})
	covs <- aperm(covs, perm = c(2,3,1))

	dimnames(covs) <- NULL
	covs
}