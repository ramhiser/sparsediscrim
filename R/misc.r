# Returns a list of mutually exclusive folds utilizing leave-k-out crossvalidation.
cv_partition <- function(y, k = 5, seed = NULL) {
	if(!is.null(seed)) {
		set.seed(seed)
	}
	n <- length(y)
	folds <- split(sample(seq_len(n), n), gl(n = ceiling(n / k), k = k, length = n))
	names(folds) <- paste("Fold", names(folds), sep = "")
	folds
}