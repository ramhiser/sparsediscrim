#' Computes the eigenvalue decomposition of the maximum likelihood estimator
#' (MLE) of the covariance matrix for the data matrix \code{x} using the Fast
#' Singular Value Decomposition (SVD).
#'
#' First, we use the \code{fast.svd} function from the \code{corpcor} package to
#' quickly compute the SVD of the matrix \code{x}. Then, we calculate the
#' eigenvalues and their corresponding eigenvectors from the SVD.
#'
#' This function employs a well-known trick for tall data (large \code{n}, small
#' \code{p}) and wide data (large \code{p}, small \code{n}) to compute the
#' eigenvalue decomposition of the covariance matrix estimator from the data
#' matrix, \code{x}.
#'
#' TODO: Briefly describe the Fast SVD trick.
#' TODO: Describe how we use the SVD to compute the eigenvalue decomposition.
#' 
#' @export
#' @param x data matrix with \code{n} observations and \code{p} feature vectors
#' @param tol tolerance value below which the singular values of \code{x} are
#' considered zero.
#' @return a list containing:
#'   the eigenvectors of the covariance matrix MLE
#'   the eigenvalues of the covariance matrix MLE
fast_cov_eigen <- function(x, tol = 1e-6) {
  require('corpcor')
  # First, we center 'x'
  x <- scale(x, center = TRUE, scale = FALSE)
  n <- nrow(x)
  svd_out <- fast.svd(x, tol = tol)
  list(vectors = svd_out$v, values = svd_out$d^2 / n)
}

#' Computes the eigenvalue decomposition of the pooled sample covariance matrix
#' for the data matrix \code{x} using the Fast Singular Value Decomposition (SVD).
#'
#' First, we use the \code{fast.svd} function from the \code{corpcor} package to
#' quickly compute the SVD of the matrix \code{x}. Then, we calculate the
#' eigenvalues and their corresponding eigenvectors from the SVD.
#'
#' This function employs a well-known trick for tall data (large \code{n}, small
#' \code{p}) and wide data (large \code{p}, small \code{n}) to compute the
#' eigenvalue decomposition of the covariance matrix estimator from the data
#' matrix, \code{x}.
#'
#' TODO: Briefly describe the Fast SVD trick.
#' TODO: Describe how we use the SVD to compute the eigenvalue decomposition.
#' 
#' @export
#' @param x data matrix with \code{n} observations and \code{p} feature vectors
#' @param y class labels for observations (rows) in \code{x}
#' @param tol tolerance value below which the singular values of \code{x} are
#' considered zero.
#' @return a list containing:
#'   the eigenvectors of the covariance matrix MLE
#'   the eigenvalues of the covariance matrix MLE
fast_cov_pool_eigen <- function(x, y, tol = 1e-6) {
  require('corpcor')
  # First, we center each class.
  x <- tapply(seq_along(y), y, function(i) {
    scale(x[i, ], center = TRUE, scale = FALSE)
  })
  x <- do.call(rbind, x)
  n <- nrow(x)
  svd_out <- fast.svd(x, tol = tol)
  list(vectors = svd_out$v, values = svd_out$d^2 / n)
}
