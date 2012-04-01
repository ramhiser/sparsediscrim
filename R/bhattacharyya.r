#' Bhattacharyya Distance for Simultaneously Diagonalized Data
#'
#' For a data matrix, \code{x}, that has been simultaneously diagonalized, we
#' calculate the Bhattacharyya distance (divergence) between the two classes
#' given in \code{y}.
#'
#' We assume that the matrix \code{x} has been simultaneously diagonalized. If
#' simultaneous diagonalization has not been employed, this function, by
#' default, will calculate the Bhattacharyya distance between the two classes
#' under the assumption that each class is multivariate normal with diagonal
#' covariance matrices. However, in this case any choice for \code{q}, other
#' than the default, will not be applicable.
#'
#' By default, we shrink the (diagonal) maximum likelihood estimators (MLEs) for
#' the covariance matrices with the Minimum Distance Empirical Bayes estimator
#' (MDEB) from Srivastava and Kubokawa (2007), which effectively shrinks the MLEs
#' towards an identity matrix that is scaled by the average of the nonzero
#' eigenvalues.
#'
#' The value of \code{q} determines the number of columns that will be used to
#' calculate the Bhattacharyya distance. This approach is useful in determining
#' an 'optimal' value for \code{q} over a range of candidate values.
#'
#' TODO: Provide a reference for Fukunaga's text for Bhattacharyya distance.
#' TODO: Provide a reference for Srivastava and Kubokawa (2007).
#'
#' @export
#' @param x data matrix with \code{n} observations and \code{p} feature vectors
#' @param y class labels for observations (rows) in \code{x}
#' @param q the number of columns of \code{x} to use in calculating the
#' Bhattacharyya distance. By default, we use \code{ncol(x)}.
#' @param shrink By default, (if \code{TRUE}), we shrink each covariance matrix
#' with the MDEB covariance matrix estimator. Otherwise, no shrinkage is applied.
#' @return the estimated Bhattacharyya distance between the two classes given in
#' \code{y}.
bhatta_simdiag <- function(x, y, q = NULL, shrink = TRUE) {
  x <- as.matrix(x)
  y <- as.factor(y)

  # Selects the first 'q' columns of 'x'. By default, we keep all of the columns
  # of 'x'.
  if (is.null(q)) {
    q <- ncol(x)
  }
  x <- x[, seq_len(q)]

  # Partitions the matrix 'x' into the data for each class.
  class_labels <- levels(y)
  x1 <- x[which(y == class_labels[1]), ]
  x2 <- x[which(y == class_labels[2]), ]

  # Calculates the MLEs for the population means and (diagonal) covariance
  # matrices for each class.
  xbar1 <- colMeans(x1)
  xbar2 <- colMeans(x2)
  diff_xbar <- as.vector(xbar1 - xbar2)
  
  cov1 <- cov_mle(x1, diag = TRUE)
  cov2 <- cov_mle(x2, diag = TRUE)

  # Shrink the diagonal covariance matrices, if specified.
  if (shrink) {
    n1 <- nrow(x1)
    n2 <- nrow(x2)

    cov1 <- cov1 + sum(cov1) / min(n1, q)
    cov2 <- cov1 + sum(cov2) / min(n2, q)
  }
  cov_avg <- (cov1 + cov2) / 2

  # Calculates the determinants of the diagonal covariance matrices.
  det_1 <- prod(cov1)
  det_2 <- prod(cov2)
  det_avg <- prod(cov_avg)

  # Calculates the estimated Bhattacharyya distance between the two classes.
  distance <- (1/8) * sum(diff_xbar^2 / cov_avg)
  distance <- distance + log(det_avg / sqrt(det_1 * det_2)) / 2
  distance
}
