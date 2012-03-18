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
