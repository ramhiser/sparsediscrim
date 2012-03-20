#' Transforms a matrix or vector with a simultaneous diagonalizing matrix.
#'
#' TODO
#'
#' @export
#' @param obj object of type 'simdiag'
#' @param x matrix or vector
#' @param q integer. (Optional) Dimension of subspace to which to project.
#' @return matrix
simdiag_transform <- function(obj, x, q = NULL) {
  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }

  if (is.null(q)) {
    q <- nrow(obj$Q)
  } else {
    q <- as.integer(q)
    if (q > nrow(obj$Q)) {
      stop("The reduced dimension 'q' cannot be larger than the number of columns of 'x'.")
    }
  }
  tcrossprod(x, obj$Q[seq_len(q), ])
}
