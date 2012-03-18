#' Simultaneously diagonalizes the sample covariance matrices for a data set
#' generated from 2 classes with the Fast SVD method.
#'
#' This function determines the matrix, \code{Q}, that simultaneously diagonalizes
#' the sample covariance matrices for a data set generated from 2 classes.
#' If the matrices are singular, we employ a dimension reduction
#' method that simultaneously diagonalizes them to a lower dimension, \code{q}.

#' This function simultaneously diagonalizes two real symmetric matrices.
#' If the matrices are singular, we employ a dimension reduction
#' method that simultaneously diagonalizes them to a lower dimension, \code{q}.
#'
#' The user can manually provides the dimension, \code{q}. This is potentially
#' useful if a scree plot is utilized or a manual threshold is desired.
#' Otherwise, we select \code{q} automatically to be the number of non-zero
#' eigenvalues. To determine what constitutes 'non-zero', we select all of the
#' eigenvalues larger than the specified tolerance, \code{tol}. This is similar
#' to a typical automatic approach to Principal Components Analysis (PCA).
#' 
#' The returned simultaneous diagonalizing matrix, \code{Q}, is constructed such
#' that
#'
#' \deqn{Q \\Sigma_1 Q' = I_q,}
#' and
#' \deqn{Q \\Sigma_2 Q' = D},
#'
#' where \eqn{I_q} denotes the \eqn{q}-dimensional identity matrix and \eqn{D}
#' is the diagonal matrix with diagonal entries as the characteristic roots
#' (generalized eigenvalues) of det(\\Sigma_1 - \\lambda \\Sigma_2) = 0.
#'
#' We use the 'Fast Singular Value Decomposition (SVD)' approach to determine the
#' eigenvectors and eigenvalues of the class with the larger rank. This technique
#' provides significant improvements in computational runtime when the number of
#' observations number is much less than the number of features. That is, we gain
#' significantly in terms of runtime if the number of rows of \code{x} is much
#' less than the number of columns of \code{x}.
#'
#' We remark that the rank estimation (i.e. choosing the value of \code{q} when
#' it is not manually specified) can be subject to numerical instabilities if the
#' value for \code{tol} is too small. In fact, the ranks of the covariance
#' matrices can be overestimated, and therefore, the value for \code{q} can be
#' too larger; this issue yields a linear transformation, \code{Q}, that projects
#' the data to a dimension larger than the actual rank. The result is that sample
#' covariance matrices might not be fully simultaneously diagonalized. Our
#' approach to overcome the issue is to set \code{q} to \code{max(n1, n2) - 1} if
#' \code{q >= max(n1, n2)} because the ranks of the class sample covariance
#' matrices are bounded above by the class sample size minus 1.
#'
#' @param x matrix containing the training data. The rows are the sample
#' observations, and the columns are the features.
#' @param y vector of class labels for each training observation
#' @param calc_ranks logical value. We consider the first class to be the class
#' with sample covariance of larger rank. By default, if \code{calc_ranks} is
#' \code{FALSE}, we consider the first class to be the one with a larger sample
#' size.in terms of the larger. Otherwise, if \code{TRUE}, we calculate the ranks
#' numerically as the number of nonzero eigenvalues and select the first class to
#' have larger rank. Note that this option adds additional computations. If
#' multicollinearity is suspected and if the sample sizes are approximately
#' equal, we recommend that \code{calc_ranks} be set to \code{TRUE}.
#' @param q the reduced dimension of the simultaneous diagonalizer. By default,
#' \code{q} is selected to be the rank of the sample covariance matrix for the
#' first class.
#' @param tol a value indicating the magnitude below which eigenvalues are
#' considered 0.
#' @return matrix of size \eqn{q \times p} that simultaneously diagonalizes the
#' sample covariance matrices (MLEs) for each class
simdiag_fastsvd <- function(x, y, calc_ranks = FALSE, q = NULL,
                            tol = 1e-6) {
  x <- as.matrix(x)
  y <- as.factor(y)
  levels_y <- levels(y)

  p <- ncol(x)
  # If a value for 'q' is provided, we check that it is between 1 and p,
  # inclusively. If not, an error is thrown.
  if (!is.null(q)) {
    q <- as.integer(q)
    if (q < 1 || q > p) {
      stop("The value for 'q' must be between 1 and p, inclusively.")
    }
  }

  # We only simultaneously diagonalize the covariance matrices for two classes.
  if (length(levels_y) < 2) {
    stop("Only 1 class is given. Use PCA instead.")
  } else if (length(levels_y) > 2) {
    stop("Currently, we do not simultaneously diagonalize more than 2 classes.")
  }

  x1 <- x[which(y == levels_y[1]), ]
  x2 <- x[which(y == levels_y[2]), ]
  n1 <- nrow(x1)
  n2 <- nrow(x2)

  # Centers both 'x1' and 'x2'. We remove the 'attr' added via 'scale.'
  x1_centered <- scale(x1, scale = FALSE)
  attr(x1_centered, 'scaled:center') <- NULL
  x2_centered <- scale(x2, scale = FALSE)
  attr(x2_centered, 'scaled:center') <- NULL

  # If 'calc_ranks' is selected, we compute the Fast SVD for both classes and
  # determine the class that has the sample covariance matrix with higher rank.
  # Otherwise, (by default) we use the class with the larger sample size; a
  # larger sample size can yield a larger rank if multi-collinearity is not an
  # issue.
  if (calc_ranks) {
    fast_svd_1 <- corpcor:::fast.svd(x1_centered / sqrt(n1), tol = tol)
    fast_svd_2 <- corpcor:::fast.svd(x2_centered / sqrt(n2), tol = tol)

    if(length(fast_svd_1$d) >= length(fast_svd_2$d)) {
      fast_svd <- fast_svd_1
    } else {
      fast_svd <- fast_svd_2
      x2 <- x1
    }
    rm(fast_svd_1)
    rm(fast_svd_2)
  } else {
    if (n1 >= n2) {
      fast_svd <- corpcor:::fast.svd(x1_centered / sqrt(n1), tol = tol)
    } else {
      fast_svd <- corpcor:::fast.svd(x2_centered / sqrt(n2), tol = tol)
      x2 <- x1
    }
  }

  # We remove the centered data matrices to free up memory; this is necessary
  # especially for "big data."
  rm(x1_centered)
  rm(x2_centered)

  Q_A <- tcrossprod(diag(1/fast_svd$d), fast_svd$v)

  # By default, (if 'q' is NULL), the value of 'q' is chosen to be the estimated
  # rank of the sample covariance matrix of 'x1'.
  if(is.null(q)) {
    q <- length(fast_svd$d)

    # The larger rank of the covariance matrices is bounded above by the larger
    # sample size minus 1. So, we perform a check for this; if 'q' turns out to
    # be at least the larger sample size, we manually set it to be equal to
    # the larger sample size minus 1.
    if (q >= max(n1, n2)) {
      q <- max(n1, n2) - 1
    }
  }

  eigen_2 <- eigen(Q_A %*% tcrossprod(cov_mle(x2), Q_A), symmetric = TRUE)

  u_2 <- eigen_2$vectors[, seq_len(q)]
  Q <- crossprod(u_2, Q_A)

  list(Q = Q, q = q)
}
