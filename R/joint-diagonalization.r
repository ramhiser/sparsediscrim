library(jointDiag)

joint.diagonalization <- function(df, method = c("none", "general-eigen", "asfari", "jade", "jedi", "qdiag", "ffdiag", "jadiag", "uwedge"),
	tol = 1e-6, max.iter = 250, shrink = TRUE, shrink.val = 0.01) {

	method <- match.arg(method)
	if(method == "none") {
		B <- diag(ncol(df) - 1)
	} else if(method == "general-eigen") {
		general.eigen.cov <- geneigen.cov(df, shrink = FALSE, shrink.val = 0.01)
		B <- diagonalize.geneigen(general.eigen.cov[[1]], general.eigen.cov[[2]])
		df <- joint.diagonalization.transform(df, B)
	} else if(method == "asfari") {
		asfari.cov <- asfari.cov(df, shrink = TRUE, shrink.val = shrink.val)
		B <- LUJID(asfari.cov, mode = 'B', ERR = tol, RBALANCE = 3, ITER = max.iter)
		df <- joint.diagonalization.transform(df, B)
	} else if(method == "jade") {
		warning("JADE algorithm has not been implemented yet.")
		B <- diag(ncol(df) - 1)
	} else if(method == "jedi") {
		ajd.cov <- jointDiag.cov(df, shrink = TRUE, shrink.val = shrink.val)
		B <- jedi(ajd.cov, eps = tol, itermax = max.iter)$A
		df <- joint.diagonalization.transform(df, B)
	} else {
		ajd.cov <- jointDiag.cov(df, shrink = TRUE, shrink.val = shrink.val)
		B <- ajd(ajd.cov, method = method, eps = tol, itermax = max.iter)$B
		df <- joint.diagonalization.transform(df, B)
	}
	list(transformed.df = df, B = B, method = method)
}

# This function computes a linear transformation of the individual matrices
# within the dataframe "df" with the matrix B.
# We assume that the first column of the data.frame df is named "labels"
# and contains the class labels.
# Returns a data.frame, where the data matrix for each class has been
# transformed by post-multiply by t(B).
joint.diagonalization.transform <- function(df, B) {
	transformed.df <- ddply(df, .(labels), function(class.df) {
		x <- as.matrix(class.df[,-1])
		dimnames(x) <- NULL
		data.frame(x %*% t(B))
	})
	transformed.df
}