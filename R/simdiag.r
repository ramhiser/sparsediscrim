#' Simultaneously diagonalizes two real symmetric matrices.
#'
#' This function simultaneously diagaonlizes two real symmetric matrices.
#' If the matrices are singular, we employ a dimension reduction
#' method that simultaneously diagonalizes them in a lower dimension.
#'
#' If the user selects the dimension reduction option, the user may
#' specify a value of k. This is potentially useful if a scree plot
#' is utilized or a manual threshold is desired. Otherwise, we
#' select k automatically to be the number of non-zero eigenvalues.
#' To determine what constitutes 'non-zero', we select all of the eigenvalues
#' larger than the specified tolerance, 'tol'.
#'
#' If the user does not wish to use the dimension reduction method,
#' then k is automatically computed as the number of eigenvalues of B.
#' Hence, in this case, k = p.
#'
#' The returned simultaneous diagonalizing matrix, Q, is constructed such that
#'
#' Q' A Q = D, and
#' Q' B Q = I_k,
#'
#' where I_k denotes the k-dimensional identity matrix and Lambda is the 
#' diagonal matrix with diagonal entries as the characteristic roots
#' (generalized eigenvalues) of det(B - lambda A) = 0.
#'
#' @param A p x p real symmetric matrix
#' @param B p x p real symmetric matrix
#' @param dim_reduce logical flag: Should we reduce the dimension of the simultaneous diagonalization?
#' @param k the reduced dimension of the simultaneous diagonalizer.
#' @param tol a value indicating the magnitude below which eigenvalues are considered 0.
#' @return a list containing:
#'	the matrix Q that is the p x k simultaneous diagonalizing matrix
#'	the selected value of k
simdiag <- function(A, B, dim_reduce = F, k = NULL, tol = sqrt(.Machine$double.eps)) {
	B_eigen <- eigen(B, symmetric = TRUE)
	
	# We note that B may singular here and denote its rank
	# as q, where 1 <= q <= p.
	# We apply our dimension reduction method that uses the top
	# k eigenvalues, where 1 <= k <= q. One approach to determining
	# k is to set it equal to the number of nonzero eigenvalues.
	# We apply the square root to R's built-in non-zero threshold
	# to declare an eigenvalue as 0.
	if(!dim_reduce) {
		k <- length(B_eigen$values)
	} else if(is.null(k)) {
		k <- sum(B_eigen$values > tol)
	}
	Q_B <- B_eigen$vectors[, seq_len(k)] %*% diag(B_eigen$values[seq_len(k)]^(-1/2))
	Q_A_eigen <- eigen(t(Q_B) %*% A %*% Q_B, symmetric = TRUE)
	list(Q = Q_B %*% Q_A_eigen$vectors, k = k)
}