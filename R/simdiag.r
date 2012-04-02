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
#' @param q the reduced dimension of the simultaneous diagonalizer.
#' @param method the method to automatically select the value of \code{q}.
#' Ignored if \code{q} is given.
#' @param simdiag_data Should the data matrix, \code{x}, be transformed and
#' returned? If TRUE (default), we return an updated matrix, \code{x}. For large
#' data sets, this option can yield a very large \code{simdiag} object and should
#' possibly be calculated manually.
#' @param calc_ranks logical value. We consider the first class to be the class
#' with sample covariance of larger rank. By default, if \code{calc_ranks} is
#' \code{FALSE}, we consider the first class to be the one with a larger sample
#' size.in terms of the larger. Otherwise, if \code{TRUE}, we calculate the ranks
#' numerically as the number of nonzero eigenvalues and select the first class to
#' have larger rank. Note that this option adds additional computations. If
#' multicollinearity is suspected and if the sample sizes are approximately
#' equal, we recommend that \code{calc_ranks} be set to \code{TRUE}. Only used if
#' the \code{method} is selected to be \code{fast_svd}; this is the default
#' option.
#' @param tol a value indicating the magnitude below which eigenvalues are
#' considered 0.
#' @return a list containing:
#' \itemize{
#'   \item \code{Q}: simultaneous diagonalizing matrix of size \eqn{p \times q}
#'	 \item \code{q}: the reduced dimension of the data
#'   \item \code{x}: the transformed data matrix, \code{x}. Returned, if
#'   \code{simdiag_data} is TRUE. Otherwise, \code{x} is not returned.
#' }
simdiag <- function(x, ...)
  UseMethod("simdiag")

#' @rdname simdiag
#' @method simdiag default
#' @S3method simdiag default
simdiag.default <- function(x, y, q = NULL,
                            fast_svd = TRUE,
                            bhattacharyya = TRUE,
                            tol = 1e-5) {
  x <- as.matrix(x)
  y <- as.factor(y)

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
      bhatta_out <- dimred_bhatta_simdiag(obj$x, y)
      obj$q <- bhatta_out$q
      obj$bhattacharyya <- bhatta_out$bhattacharyya
      obj$pct_change <- bhatta_out$pct_change
    }
  }

  obj$Q <- obj$Q[seq_len(obj$q), ]
  obj$x <- obj$x[, seq_len(obj$q)]               

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
simdiag.formula <- function(formula, data, q = NULL,
                            method = c("rank", "bhattacharyya"), tol = 1e-6) {
  # The formula interface includes an intercept. If the user includes the
  # intercept in the model, it should be removed. Otherwise, errors and doom
  # happen.
  # To remove the intercept, we update the formula, like so:
  # (NOTE: The terms must be collected in case the dot (.) notation is used)
  formula <- diagdiscrim:::no_intercept(formula, data)
  
  mf <- model.frame(formula = formula, data = data)
  x <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)

  obj <- simdiag.default(x, y, q = q, method = method, tol = tol)
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

