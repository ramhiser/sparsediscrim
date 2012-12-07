#' Transforms a matrix or vector with a simultaneous diagonalizing matrix.
#'
#' TODO
#'
#' @export
#' @param obj object of type 'simdiag'
#' @param x matrix or vector
#' @return matrix
simdiag_transform <- function(object, x) {
  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }

  tcrossprod(x, object$Q)
}
