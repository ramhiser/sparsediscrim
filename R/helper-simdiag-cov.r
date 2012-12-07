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
#' @export
#' @param A p x p real symmetric matrix
#' @param B p x p real symmetric matrix
#' @param q the reduced dimension of the simultaneous diagonalizer.
#' @param tol a value indicating the magnitude below which eigenvalues are
#' considered 0.
#' @return a list containing:
#'	the matrix Q that is the \eqn{p \times q} simultaneous diagonalizing matrix
#'	the selected value of \code{q}
simdiag_cov <- function(A, B, q = NULL, tol = 1e-6) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  
  if (nrow(A) != ncol(A)) {
    stop("The matrix, 'A', must be square.")
  }

  if (nrow(B) != ncol(B)) {
    stop("The matrix, 'B', must be square.")
  }
  if (nrow(A) != nrow(B)) {
    stop("The matrices, 'A' and 'B', must have the same size.")
  }

  # If a value for 'q' is provided, we check that it is between 1 and p,
  # inclusively. If not, an error is thrown.
  p <- nrow(A)
  if (!is.null(q)) {
    q <- as.integer(q)
    if (q < 1 || q > p) {
      stop("The value for 'q' must be between 1 and p, inclusively.")
    }
  }
  
  eigen_A <- eigen(A, symmetric = TRUE)
  if (is.null(q)) {
    q <- sum(eigen_A$values > tol)
  }
  U1 <- eigen_A$vectors[, seq_len(q)]
  Lambda_1_inv <- diag(eigen_A$values[seq_len(q)]^(-1/2))
  Q1 <- tcrossprod(lambda_1_inv, U1)

  eigen_2 <- eigen(Q1 %*% tcrossprod(B, Q1), symmetric = TRUE)
  U2 <- eigen_2$vectors[, seq_len(q)]
  Q <- crossprod(U2, Q1)

  list(Q = Q, q = q)
}
