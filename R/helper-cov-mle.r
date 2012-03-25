#' Computes the maximum likelihood estimator for the sample covariance matrix
#' under the assumption of multivariate normality.
#'
#' For a sample matrix, \code{x}, we compute the sample covariance matrix of the
#' data as the maximum likelihood estimator (MLE) of the population covariance
#' matrix.
#'
#' @export
#' @param x data matrix with \code{n} observations and \code{p} feature vectors
#' @return sample covariance matrix of size \eqn{p \times p}
cov_mle <- function(x) {
  (nrow(x) - 1) / nrow(x) * cov(x)
}

#' Computes the pooled maximum likelihood estimator for the common covariance
#' matrix from \eqn{K} multivariate normal populations.
#'
#' TODO: Description
#'
#' @export
#' @param x data matrix with \code{n} observations and \code{p} feature vectors
#' @param y class labels for observations (rows) in \code{x}
#' @return pooled sample covariance matrix of size \eqn{p \times p}
#' @examples
#' cov_pool(iris[,-5], iris$Species)
cov_pool <- function(x, y) {
  x <- as.matrix(x)
  y <- as.factor(y)
  n <- length(y)
  scatter_matrices <- tapply(seq_len(n), y, function(i) {
    (length(i) - 1) * cov(x[i, ])
  })
  Reduce("+", scatter_matrices) / n
}
