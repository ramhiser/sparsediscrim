#' Simultaneously diagonalizes the sample covariance matrices for a data set
#' generated from \code{K} classes.
#'
#' This function determines the matrix, \code{Q}, that simultaneously
#' diagonalizes the sample covariance matrices for a data set generated from
#' \code{K} classes. If the matrices are singular, we employ a dimension
#' reduction method that simultaneously diagonalizes them to a lower dimension,
#' \code{q}.
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
#' \deqn{Q \widehat{\\Sigma_k} Q' = D_k}
#'
#' for all \eqn{k = 1, \ldots, K}, where \eqn{\widehat{\\Sigma_k}} is the maximum
#' likelihood estimator (under normality) for the \eqn{k}th class covariance
#' matrix, \eqn{\\Sigma_k}, where \eqn{D_k} is a diagonal matrix.
#'
#' The returned simultaneously-diagonalizing matrix, \code{Q}, is of size
#' \eqn{p \times q}.
#'
#' Currently, we only allow \eqn{K = 2}. The \eqn{K = 1} case is equivalent to
#' PCA. The \eqn{K \ge 2} case is of ongoing research; note that a set of
#' \code{K} positive semidefinite matrices are simultaneously diagonalizable if
#' they are pairwise-commute (i.e. \eqn{A B = B A} for any two matrices \eqn{A}
#' and \eqn{B} in the set). However, pairwise commuting is a strict assumption
#' that we wish to relax.
#'
#' @export
#' @param x matrix containing the training data. The rows are the sample
#' observations, and the columns are the features.
#' @param y vector of class labels for each training observation
#' @param q the reduced dimension of the simultaneous diagonalizer. If
#' \code{NULL} (default), we choose \code{q} via the Bhattacharyya distance
#' method, discussed below.
#' @param fast_svd logical value that indicates if we should utilize the Fast SVD
#' compute the eigenvalue decomposition of the sample covariance matrix having
#' larger rank. By default, we use this Fast SVD method, which is computationally
#' far superior for 'wide data' (large \code{p}, small {n}). 
#' @param bhattacharyya logical value that indicates if we should select the
#' reduced dimension, \code{q}, via our proposed Bhattacharyya method. If
#' \code{TRUE}, we choose the value for \code{q} to be the dimension, where the
#' Bhattacharyya distance between the two classes hardly changes for
#' \code{q + 1}. Our technique is similar, in spirit, to the 'elbow criterion' in
#' a Principal Components Analysis scree plot. The \code{bhattacharyya} argument
#' is ignored if a value for \code{q} is specified.
#' @param tol a value indicating the magnitude below which eigenvalues are
#' considered 0.
#' @param bhatta_prop threshold value for the maximum cumulative proportion of
#' the Bhattacharyya distance. We use the threshold to determine the reduced
#' dimension \code{q}. When paired with the SimDiag method, the Bhattacharyya
#' approach to dimension reduction generalizes the scree plot idea that is often
#' utilized with Principal Components Analysis. Ignored if \code{bhattacharyya}
#' is \code{FALSE}.
#' @param bhatta_shrink By default, (if \code{TRUE}), we shrink each covariance
#' matrix with the MDEB covariance matrix estimator. Otherwise, no shrinkage is
#' applied.  Ignored if \code{bhattacharyya} is \code{FALSE}.

#' @return a list containing:
#' \itemize{
#'   \item \code{Q}: simultaneous diagonalizing matrix of size \eqn{p \times q}
#'	 \item \code{q}: the reduced dimension of the data
#'   \item \code{x}: the transformed data matrix, \code{x}.
#' }
simdiag <- function(x, ...)
  UseMethod("simdiag")

#' @rdname simdiag
#' @method simdiag default
#' @S3method simdiag default
simdiag.default <- function(x, y, q = NULL, bhattacharyya = TRUE,
                            fast_svd = TRUE, tol = 1e-6, bhatta_prop = 0.9,
                            bhatta_shrink = FALSE) {
  x <- as.matrix(x)
  y <- as.factor(y)
  p <- ncol(x)

  # For now, we only simultaneously diagonalize the covariance matrices for two
  # classes.
  if (nlevels(y) < 2) {
    stop("Only 1 class is given. Use PCA instead.")
  } else if (nlevels(y) > 2) {
    stop("Currently, we do not simultaneously diagonalize more than 2 classes.")
  }

  # If a value for 'q' is provided, we check that it is between 1 and p,
  # inclusively. If not, an error is thrown.
  if (!is.null(q)) {
    q <- as.integer(q)
    if (q < 1 || q > p) {
      stop("The value for 'q' must be between 1 and p, inclusively.")
    }
  }

  # Computes the MLEs of the covariance matrices for each class
  # Because we assume rank(Sig1) >= rank(Sig2), we designate class 1 as the class
  # with the larger sample size, which is our preliminary candidate value for 'q'.
  class_label1 <- levels(y)[which.max(table(y))]
  class_label2 <- levels(y)[which.min(table(y))]

  x1 <- x[which(y == class_label1), ]
  x2 <- x[which(y == class_label2), ]

  cov1 <- cov_mle(x1)
  cov2 <- cov_mle(x2)

  if (fast_svd) {
    cov1_eigen <- fast_cov_eigen(x1)
  } else {
    cov1_eigen <- eigen(cov1, symmetric = TRUE)
  }
  
  Lambda1 <- cov1_eigen$values
  q1 <- sum(Lambda1 >= tol)
  U1 <- cov1_eigen$vectors[, seq_len(q1)]
  Lambda1_pseudo <- 1 / Lambda1[seq_len(q1)]
  Lambda1_pseudo_sqrt <- Lambda1_pseudo^(1/2)
  Q1 <- tcrossprod(diag(Lambda1_pseudo_sqrt), U1)

  eigen_out <- eigen(Q1 %*% tcrossprod(cov2, Q1), symmetric = TRUE)
  U2 <- eigen_out$vectors
  Lambda2 <- eigen_out$values

  # Creates a list named 'obj' that will contain all of the returned information.
  obj <- list()

  obj$Q <- crossprod(U2, Q1)

  # We transform the matrix 'x' and store it in the returned named list as 'x'.
  obj$x <- tcrossprod(x, obj$Q)

  # By default, we select the value of 'q' to be the number of rows of Q. If the
  # user has specified a different value for 'q', we use it unless it exceeds
  # the number of rows of Q, in which case we use the default value.
  obj$q <- nrow(obj$Q)
  if (!is.null(q)) {
    if (q <= nrow(obj$Q)) {
      obj$q <- q
    } else {
      warning("The value for 'q' exceeds the number of rows of 'Q'. Using 'nrow(Q)' instead.")
    }
  } else {
    if (bhattacharyya) {
      bhatta_out <- bhatta_simdiag(x, y, bhatta_prop = bhatta_prop,
                                   shrink = bhatta_shrink)
      obj$q <- bhatta_out$q
      obj$bhattacharyya <- bhatta_out$bhattacharyya
      obj$pct_change <- bhatta_out$pct_change
    }
  }

  obj$Q <- obj$Q[seq_len(obj$q), ]
  obj$x <- obj$x[, seq_len(obj$q)]
  obj$Lambda1 <- Lambda1
  obj$Lambda2 <- Lambda2
  
  # Creates an object of type 'simdiag' and adds the 'match.call' to the object
  obj$call <- match.call()
  class(obj) <- "simdiag"
	
	obj
}

#' @param formula formula of the form \code{groups ~ x1 + x2 + ...} That is,
#' the response is the grouping factor and the right hand side specifies the
#' (non-factor) discriminators.
#' @param data data frame from which variables specified in \code{formula} are
#' preferentially to be taken.
#' @rdname simdiag
#' @method simdiag formula
#' @S3method simdiag formula
simdiag.formula <- function(formula, data, ...) {
  # The formula interface includes an intercept. If the user includes the
  # intercept in the model, it should be removed. Otherwise, errors and doom
  # happen.
  # To remove the intercept, we update the formula, like so:
  # (NOTE: The terms must be collected in case the dot (.) notation is used)
  formula <- no_intercept(formula, data)
  
  mf <- model.frame(formula = formula, data = data)
  x <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)

  obj <- simdiag.default(x, y, ...)
  obj$call <- match.call()
  obj$formula <- formula
  obj
}

#' Outputs the summary for a 'simdiag' object.
#'
#' Summarizes the simultaneous diagonalization information for a set of classes.
#'
#' @keywords internal
#' @param x object to print
#' @param ... unused
#' @rdname simdiag
#' @method print simdiag
#' @S3method print simdiag
#' @export
print.simdiag <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
}

