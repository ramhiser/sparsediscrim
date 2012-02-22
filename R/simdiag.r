#' Simultaneously diagonalizes the sample covariance matrices for a data set
#' generated from \code{K} classes.
#'
#' This function determines the matrix, \code{Q}, that simultaneously diagonalizes
#' the sample covariance matrices for a data set generated from \code{K} classes.
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
#' \deqn{Q \widehat{\Sigma_k} Q' = D_k}
#'
#' for all \eqn{k = 1, \ldots, K}, where \eqn{\widehat{\Sigma_k}} is the maximum
#' likelihood estimator (under normality) for the \eqn{k}th class covariance
#' matrix, \eqn{\Sigma_k}, where \eqn{D_k} is a diagonal matrix.
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
simdiag.default <- function(x, y, q = NULL, method = c("rank", "bhattacharyya"),
                            simdiag_data = TRUE, tol = 1e-6) {
  x <- as.matrix(x)
  y <- as.factor(y)

  p <- ncol(x)
  # If a value for 'q' is provided, we check that it is between 1 and p,
  # inclusively. If not, an error is thrown.
  if (!is.null(q)) {
    q <- as.integer(q)
    if (q < 1 || q > p) {
      stop("The value for 'q' must be between 1 and p, inclusively.")
    }
  }

  # For now, we only simultaneously diagonalize the covariance matrices for two
  # classes.
  if (nlevels(y) < 2) {
    stop("Only 1 class is given. Use PCA instead.")
  } else if (nlevels(y) > 2) {
    stop("Currently, we cannot simultaneously diagonalize more than 2 classes.")
  }

  # Computes the MLEs of the covariance matrices for each class
  # Because we assume rank(Sig1) >= rank(Sig2), we designate class 1 as the class
  # with the larger sample size, which is our preliminary candidate value for 'q'.
  larger_class <- which.max(table(y))
  n1 <- sum(y == larger_class)
  n2 <- sum(y != larger_class)
  Sig1 <- (n1 - 1) * cov(x[which(y == larger_class), ]) / n1
  Sig2 <- (n2 - 1) * cov(x[which(y != larger_class), ]) / n2

  # In future work, we intend to allow for more than the K = 2 case. Hence, we
  # have delegated the actual simultaneous diagonalization to the function,
  # 'simdiag2'.
  # TODO: Pass 'method' to simdiag2. For now, we are only using the rank
  # estimation to select 'q'.
  obj <- diagdiscrim:::simdiag2(A = Sig1, B = Sig2, q = q, tol = tol)

  # Creates an object of type 'simdiag' and adds the 'match.call' to the object
  obj$call <- match.call()
  class(obj) <- "simdiag"

  # If the user has specified to simultaneously diagonalize the data, we
  # transform the matrix 'x' and store it in the returned named list as 'x'.
  # NOTE: This is not recommended if limited memory is an issue, or if the data
  # set is extremely large.
  if (simdiag_data) {
    obj$x <- simdiag_transform(obj, x)
  }
	
	obj
}

#' @param formula A formula of the form \code{groups ~ x1 + x2 + ...} That is,
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

