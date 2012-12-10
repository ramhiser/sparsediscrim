#' Simultaneously diagonalizes two real symmetric matrices.
#'
#' This function simultaneously diagonalizes two real symmetric matrices.
#' If the matrices are singular, we employ a dimension reduction
#' method that simultaneously diagonalizes them to a lower dimension, \code{q}.
#'
#' The user can manually provides the dimension, \code{q}. This is potentially
#' useful if a scree plot is utilized or a manual threshold is desired.
#' Otherwise, we select \code{q} automatically to be the number of non-zero
#' eigenvalues. To determine what constitutes 'non-zero', we select all of the
#' eigenvalues larger than the specified tolerance, \code{tol}. This is similar
#' to a typical automatic approach to Principal Components Analysis.
#' 
#' The returned simultaneous diagonalizing matrix, \code{Q}, is constructed such
#' that
#'
#' \deqn{Q A Q' = D,}
#' and
#' \deqn{Q B Q' = I_q},
#'
#' where \eqn{I_q} denotes the \eqn{q}-dimensional identity matrix and \eqn{D} is
#' the diagonal matrix with diagonal entries as the characteristic roots
#' (generalized eigenvalues) of det(B - \\lambda A) = 0.
#'
#' TODO: Discuss briefly that we set 'q' to the number of nonzero eigenvalues if
#' the Fast SVD is used.
#'
#' @export
#' @param x matrix containing the training data. The rows are the sample
#' observations, and the columns are the features.
#' @param y vector of class labels for each training observation
#' @param fast_svd TODO
#' @param q the reduced dimension of the simultaneous diagonalizer.
#' @param reduced_rank Should we select \code{q} based on the rank of the
#' sample covariance matrices? Ignored if \code{q} is given. See details.
#' @param shrink TODO
#' @param tol a value indicating the magnitude below which eigenvalues are
#' considered 0.
#' @return a list containing:
#'	the matrix Q that is the \eqn{p \times q} simultaneous diagonalizing matrix
#'  the matrix x consisting of the simultaneously diagonalized data
#'	the selected value of \code{q}
simdiag_cov <- function(x, y, fast_svd = FALSE, q = NULL, reduced_rank = TRUE,
                        shrink = FALSE, tol = 1e-6) {
  x <- as.matrix(x)
  y <- as.factor(y)
  p <- ncol(x)

  # For now, we only simultaneously diagonalize the covariance matrices for two
  # classes.
  if (nlevels(y) < 2) {
    stop("Only 1 class is given. Use PCA instead.")
  } else if (nlevels(y) > 2) {
    stop("Currently, we do not simultaneously diagonalize more than 2 classes.")
  }

  # Computes the sample covariance matrices of each class.
  class_larger <- levels(y)[which.max(table(y))]
  x1 <- x[y == class_larger, ]
  x2 <- x[y != class_larger, ]

  # If the Fast SVD option is selected, we compute the eigenvalue decomposition
  # of the first class' sample covariance matrix. Otherwise, we compute the
  # eigenvalue decomposition manually.  
  if (fast_svd) {
    cov1_eigen <- fast_cov_eigen(x1, tol = tol)
    # If 'q' < the number of positive eigenvalues when the Fast SVD is used,
    # we set 'q' to be the number of positive eigenvalues.   
    if (!is.null(q) && length(cov1_eigen$values) < q) {
      warning("The reduced dimension 'q' exceeds the number of positive eigenvalues.")
      q <- length(cov1_eigen$values)
    }
  } else {
    cov1 <- cov_mle(x1)
    cov1_eigen <- eigen(cov1, symmetric = TRUE)
  }

  # Eigenvalue decomposition of the first sample covariance matrix.
  Lambda1 <- cov1_eigen$values
  U1 <- cov1_eigen$vectors

  # Now, we select the reduced dimension, 'q'. Three options:
  # 1. If a value for 'q' is provided, we check that it is between 1 and p,
  # inclusively; if not, an error is thrown.
  # 2. Select 'q' as the rank of the first sample covariance matrix. The rank is
  # calculated as the number of eigenvalues exceeding the tolerance value 'tol'
  # to compute numerically the number of nonzero eigenvalues.
  # 3. Do not apply dimension reduction.
  if (!is.null(q)) {
    q <- as.integer(q)
    if (q < 1 || q > p) {
      stop("The value for 'q' must be between 1 and p, inclusively.")
    }
  } else if (reduced_rank) {
    q <- sum(Lambda1 > tol)
  } else {
    # No dimension reduction.
    q <- p
  }

  # Reduce the rank of the eigenvalue decomposition.
  Lambda1 <- Lambda1[seq_len(q)]
  U1 <- U1[, seq_len(q)]

  # We set the inverse of the zero eigenvalues to 0 as done with the
  # Moore-Penrose pseudoinverse.
  Lambda1_pseudo_sqrt <- Lambda1^(-1/2)
  Lambda1_pseudo_sqrt <- replace(Lambda1_pseudo_sqrt, Lambda1 < tol, 0)
  
  # Whitening transform in a q-dimensional subspace
  Q1 <- tcrossprod(diag(Lambda1_pseudo_sqrt), U1)

  # Construct the simultaneous diagonalizer 'Q by diagonalizing
  # Q1 %*% Sigma2 %*% t(Q1)
  cov2 <- cov_mle(x2)
  eigen_out <- eigen(Q1 %*% tcrossprod(cov2, Q1), symmetric = TRUE)
  U2 <- eigen_out$vectors
  Lambda2 <- eigen_out$values

  # The simultaneous diagonalizer 'Q'
  Q <- crossprod(U2, Q1)

  # Simultaneously diagonalizes the data set 'x'
  x <- tcrossprod(x, Q)

  list(Q = Q, x = x, q = q)
}
