#' Tests if a matrix is diagonal.
#'
#' This function tests if the matrix \code{x} is diagonal by computing the
#' Frobenius norm of the difference between \code{x} and the matrix consisting of
#' the diagonal elements of \code{x} is less than some tolerance.
#'
#' @export
#' @param x matrix
#' @param tol tolerance value
#' @return logical value indicating whether the matrix \code{x} is diagonal
is_diagonal <- function(x, tol = 1e-8) {
  norm(diag(diag(x)) - x, "F") < tol
}
