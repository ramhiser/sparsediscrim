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

#' Determines the reduced dimension of the simultaneously diagonalized data
#' based on the Bhattacharyya distance.
#'
#' For a data matrix, \code{x}, that has been simultaneously diagonalized, we
#' calculate the Bhattacharyya distance (divergence) between the two classes
#' given in \code{y} for the dimensions ranging from 2 to \code{ncol(x)}.
#' We assume that the matrix \code{x} has been simultaneously diagonalized.
#'
#' We wish to determine the value of \code{q}, such that the percentage change
#' in the Bhattacharyya distance for the first \code{q} columns and the
#' Bhattacharyya distance for the first \code{q + 1} columns is less than
#' \code{tol}. The selected value of \code{q} is the dimension to which the data
#' matrix, \code{x}, is reduced. If are no cases where the percentage change is
#' less than \code{tol}, then we select the reduced dimension to be the
#' \code{ncol(x)}. Notice that our procedure is similar to that of the 'elbow
#' criterion' often employed with a scree plot in a Principal Components
#' Analysis (PCA). Our approach allows us to generalize the PCA dimension
#' reduction technique, accordingly.
#'
#' By default, we shrink the (diagonal) maximum likelihood estimators (MLEs) for
#' the covariance matrices with the Minimum Distance Empirical Bayes estimator
#' (MDEB) from Srivastava and Kubokawa (2007), which effectively shrinks the MLEs
#' towards an identity matrix that is scaled by the average of the nonzero
#' eigenvalues.
#'
#' TODO: Provide a reference for Fukunaga's text for Bhattacharyya distance.
#' TODO: Provide a reference for Srivastava and Kubokawa (2007).
#'
#' @export
#' @param x data matrix with \code{n} observations and \code{p} feature vectors
#' @param y class labels for observations (rows) in \code{x}
#' @param tol tolerance value. If the percentage change in the Bhattacharyya
#' distance for a given value of \code{q} and \code{q + 1} is less than
#' \code{tol}, then we select the value of \code{q} as the optimal value. If
#' there are no cases where the percentage change is less than \code{tol}, then
#' we select the reduced dimension to be the \code{ncol(x)}.
#' @param shrink By default, (if \code{TRUE}), we shrink each covariance matrix
#' with the MDEB covariance matrix estimator. Otherwise, no shrinkage is applied.
#' @return a list containing:
#'	the selected value of \code{q}
#'  the estimated Bhattacharyya distance between the two classes for each
#' candidate value of \code{q} from 2 to \code{ncol(x)}.
#'  the percentage change between the Bhattacharyya distances for the sequence of
#' candidate values of \code{q}.
dimred_bhatta_simdiag <- function(x, y, shrink = TRUE, tol = 1e-3) {
  bhatta_q <- sapply(seq.int(2, ncol(x)), function(q) {
    bhatta_simdiag(x = x, y = y, q = q, shrink = shrink)
  })
  bhatta_pct_change <- diagdiscrim:::pct_change(bhatta_q)

  small_pct_change <- which(bhatta_pct_change < tol)
  if (length(small_pct_change) == 0) {
    q <- ncol(x)
  } else {
    q <- small_pct_change[1] + 1
  }

  list(q = q, bhattacharyya = bhatta_q, pct_change = bhatta_pct_change)
}
