#' Quadratic Form of a matrix and a vector
#'
#' We compute the quadratic form of a vector and a matrix in an
#' efficient manner. Let x be a real vector of length p, and let
#' A be a p x p real matrix. Then, we compute the quadratic form
#' q by
#'
#' q = x' A x.
#'
#' A naive way to compute the quadratic form is to explicitly write
#'
#' t(x) %*% A %*% x,
#'
#' but for large p, this operation is inefficient. We provide a
#' more efficient method below.
#'
#' Note that we have adapted the code from:
#' \url{http://tolstoy.newcastle.edu.au/R/help/05/11/14989.html}
#' 
#' @param A matrix of dimension p x p
#' @param x vector of length p
#' @return quadratic form x' A x
quadform <- function(A, x) {
  drop(crossprod(x, A %*% x))
}

#' Quadratic Form of the inverse of a matrix and a vector
#'
#' We compute the quadratic form of a vector and the inverse of a matrix in an
#' efficient manner. Let x be a real vector of length p, and let
#' A be a p x p nonsingular matrix. Then, we compute the quadratic form
#' q by
#'
#' q = x' A^{-1} x.
#'
#' A naive way to compute the quadratic form is to explicitly write
#'
#' t(x) %*% solve(A) %*% x,
#'
#' but for large p, this operation is inefficient. We provide a
#' more efficient method below.
#'
#' Note that we have adapted the code from:
#' \url{http://tolstoy.newcastle.edu.au/R/help/05/11/14989.html}
#' 
#' @param A matrix that is p x p and nonsingular
#' @param x vector of length p
#' @return quadratic form x' A^{-1} x
quadform_inv <- function(A, x) {
  drop(crossprod(x, solve(A, x)))
}

#' Centers the observations in a matrix by their respective class sample means
#'
#' @param x matrix containing the training data. The rows are the sample
#' observations, and the columns are the features.
#' @param y vector of class labels for each training observation
#' @return matrix with observations centered by its corresponding class sample
#' mean
center_data <- function(x, y) {
  x <- as.matrix(x)
  y <- as.factor(y)

  # Notice that the resulting centered data are sorted by class and do not
  # preserve the original ordering of the data.
  x_centered <- tapply(seq_along(y), y, function(i) {
    scale(x[i, ], center = TRUE, scale = FALSE)
  })
  x_centered <- do.call(rbind, x_centered)

  # Sorts the centered data to preserve the original ordering.
  orig_ordering <- do.call(c, tapply(seq_along(y), y, identity))
  x_centered[orig_ordering, ] <- x_centered
  x_centered
}
