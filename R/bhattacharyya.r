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
  if (!is.null(q)) {
    q <- as.integer(q)
    if (q < 1 || q > ncol(x)) {
      stop("The value of 'q' must be between 1 and 'ncol(x)', inclusively.")
    }
    x <- as.matrix(x[, seq_len(q)])
  }

  # Partitions the matrix 'x' into the data for each class.
  class_labels <- levels(y)
  # The calls to 'as.matrix' are to maintain a matrix form if the number of
  # columns of 'x' is 1.
  x1 <- as.matrix(x[which(y == class_labels[1]), ])
  x2 <- as.matrix(x[which(y == class_labels[2]), ])

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
#' We wish to determine the \code{q} largest Bhattacharyya distances from the
#' \code{p} features in \code{x}. We assume that the class covariance matrices
#' for both classes given in \code{y} are diagonal. The selected value of
#' \code{q} is the dimension to which the data matrix, \code{x}, is reduced.
#' Notice that our procedure is similar to the reduced dimension selection
#' technique often employed with a scree plot in a Principal Components
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
#' @param bhatta_prop threshold value for the maximum cumulative proportion of
#' the Bhattacharyya distance. We use the threshold to determine the reduced
#' dimension \code{q}. When paired with the SimDiag method, the Bhattacharyya
#' approach to dimension reduction generalizes the scree plot idea that is often
#' utilized with Principal Components Analysis.
#' @param shrink By default, (if \code{TRUE}), we shrink each covariance matrix
#' with the MDEB covariance matrix estimator. Otherwise, no shrinkage is applied.
#' @return a list containing:
#' \itemize{
#'   \item \code{q}: the reduced dimension determined via \code{bhatta_pct}
#'   \item \code{bhatta_dist}: the estimated Bhattacharyya distance between the
#'   two classes for each dimension in the transformed space.
#'   \item \code{bhatta_dist_rank}: the indices of the columns of \code{x}
#'   corresponding to the \code{q} largest Bhattacharyya distances.
#'   \item \code{bhatta_cumprop}: the cumulative proportion of the sorted
#'   Bhattacharyya distances for each column given in \code{x}.
#' }
dimred_bhatta_simdiag <- function(x, y, bhatta_prop = 0.9, shrink = TRUE) {
  bhatta_dist <- sapply(seq_len(ncol(x)), function(j) {
    bhatta_simdiag(x[, j], y, shrink = shrink)
  })
  sorted_bhatta_dist <- sort(bhatta_dist, decreasing = TRUE)
  bhatta_cumprop <- cumsum(sorted_bhatta_dist) / sum(sorted_bhatta_dist)

  q <- sum(bhatta_cumprop <= bhatta_prop)
  bhatta_dist_rank <- order(bhatta_dist, decreasing = TRUE)[seq_len(q)]

  list(q = q, bhatta_dist = bhatta_dist, bhatta_dist_rank = bhatta_dist_rank,
       bhatta_cumprop = bhatta_cumprop)
}



