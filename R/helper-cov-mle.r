#' Computes the maximum likelihood estimator for the sample covariance matrix
#' under the assumption of multivariate normality.
#'
#' For a sample matrix, \code{x}, we compute the sample covariance matrix of the
#' data as the maximum likelihood estimator (MLE) of the population covariance
#' matrix.
#'
#' If the \code{diag} option is set to \code{TRUE}, then we assume the population
#' covariance matrix is diagonal, and the MLE is computed under this assumption.
#' In this case, we return a vector of length \code{p} instead.
#'
#' @export
#' @param x data matrix with \code{n} observations and \code{p} feature vectors
#' @param diag logical value. If TRUE, assumes the population covariance matrix
#' is diagonal. By default, we assume that \code{diag} is \code{FALSE}.
#' @return sample covariance matrix of size \eqn{p \times p}. If \code{diag} is
#' \code{TRUE}, then a vector of length \code{p} is returned instead.
cov_mle <- function(x, diag = FALSE) {
  n <- nrow(x)
  if (diag) {
    (n - 1) / n * apply(x, 2, var)
  } else {
    (n - 1) / n * cov(x)
  }
}

#' Computes the covariance matrix maximum likelihood estimators for each class
#' and returns a list.
#'
#' For a sample matrix, \code{x}, we compute the sample covariance matrix for
#' each class given in the vector, \code{y}.
#' @export
#' @param x data matrix with \code{n} observations and \code{p} feature vectors
#' @param y class labels for observations (rows) in \code{x}
#' @return list of the sample covariance matrices of size \eqn{p \times p} for
#' each class given in \code{y}.
cov_list <- function(x, y) {
  x <- as.matrix(x)
  y <- as.factor(y)
  tapply(seq_along(y), y, function(i) {
    cov_mle(x[i, ])
  })
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

#' Computes a shrunken version of the maximum likelihood estimator for the
#' sample covariance matrix under the assumption of multivariate normality.
#'
#' For a sample matrix, \code{x}, we compute the sample covariance matrix
#' as the maximum likelihood estimator (MLE) of the population covariance
#' matrix and shrink it towards its diagonal.
#'
#' Let \eqn{\widehat{\Sigma}} be the MLE of the covariance matrix \eqn{\Sigma}.
#' Then, we shrink the MLE towards its diagonal by computing
#' \deqn{\widehat{\Sigma}(\gamma) = \gamma \widehat{\Sigma} + (1 - \gamma) \widehat{\Sigma} \circ I_p},
#' where \eqn{\circ} denotes the Hadamard product and \eqn{\gamma \in [0,1]}.
#'
#' For \eqn{\gamma < 1}, the resulting shrunken covariance matrix estimator is
#' positive definite, and for \eqn{\gamma = 1}, we simply have the MLE, which can
#' potentially be positive semidefinite (singular).
#'
#' The estimator given here is based on Section 18.3.1 of the text,
#' "Elements of Statistical Learning."
#'
#' TODO: Cite the text.
#'
#' @export
#' @param x data matrix with \code{n} observations and \code{p} feature vectors
#' @param gamma the shrinkage parameter. Must be between 0 and 1, inclusively.
#' By default, the shrinkage parameter is 1, which simply yields the MLE.
#' @return shrunken sample covariance matrix of size \eqn{p \times p}
cov_shrink_diag <- function(x, gamma = 1) {
  if (gamma < 0 || gamma > 1) {
    stop("The value of 'gamma' must be between 0 and 1, inclusively.")
  }
  cov_x <- diagdiscrim:::cov_mle(x)
  gamma * cov_x + (1 - gamma) * diag(diag(cov_x))
}

