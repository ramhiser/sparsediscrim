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
#' The eigen_pct is the percentage of eigenvalues to use, by which we mean
#' the cumulative sum of the eigenvalues (scaled by the sum of the eigenvalues).
#' This is similar to the typical automatic approach to Principal Components Analysis.
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
#' @export
#' @param A p x p real symmetric matrix
#' @param B p x p real symmetric matrix
#' @param dim_reduce logical flag: Should we reduce the dimension of the simultaneous diagonalization?
#' @param q the reduced dimension of the simultaneous diagonalizer.
#' @param method the method to automatically select the value of \code{q}.
#' Ignored if \code{q} is given.
#' @param tol a value indicating the magnitude below which eigenvalues are considered 0.
#' @return a list containing:
#'	the matrix Q that is the p x k simultaneous diagonalizing matrix
#'	the selected value of k
simdiag2 <- function(x, y, q = NULL, tol = 1e-6) {
  p <- ncol(x)
  # If a value for 'q' is provided, we check that it is between 1 and p,
  # inclusively. If not, an error is thrown.
  if (!is.null(q)) {
    q <- as.integer(q)
    if (q < 1 || q > p) {
      stop("The value for 'q' must be between 1 and p, inclusively.")
    }
  }

  # For now, we only simultaneously diagonalize the covariance matrices for two
  # classes.
  if (nlevels(y) < 2) {
    stop("Only 1 class is given. Use PCA instead.")
  } else if (nlevels(y) > 2) {
    stop("Currently, we cannot simultaneously diagonalize more than 2 classes.")
  }

  # Computes the MLEs of the covariance matrices for each class
  # Because we assume rank(Sig1) >= rank(Sig2), we designate class 1 as the class
  # with the larger sample size, which is our preliminary candidate value for 'q'.
  larger_class <- which.max(table(y))
  n1 <- sum(y == larger_class)
  n2 <- sum(y != larger_class)
  Sig1 <- (n1 - 1) * cov(x[which(y == larger_class), ]) / n1
  Sig2 <- (n2 - 1) * cov(x[which(y != larger_class), ]) / n2

  eigen_Sig1 <- eigen(Sig1, symmetric = TRUE)
  if (is.null(q)) {
    q <- sum(eigen_Sig1$values > tol)
  }
  u1 <- eigen_Sig1$vectors[, seq_len(q)]
  lambda1_inv <- diag(eigen_Sig1$values[seq_len(q)]^(-1/2))
  Q1 <- tcrossprod(lambda1_inv, u1)
  diag(Q1 %*% tcrossprod(Sig1, Q1))

  eigen_2 <- eigen(Q1 %*% tcrossprod(Sig2, Q1), symmetric = TRUE)
  u12 <- eigen_2$vectors[, seq_len(q)]
  Q <- crossprod(u12, Q1)

  list(Q = Q, q = q)
}
