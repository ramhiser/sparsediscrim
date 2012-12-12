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
#' @param pct threshold value for the maximum cumulative proportion of
#' the Bhattacharyya distance. We use the threshold to determine the reduced
#' dimension \code{q}. When paired with the SimDiag method, the Bhattacharyya
#' approach to dimension reduction generalizes the scree plot idea that is often
#' utilized with Principal Components Analysis. Ignored if \code{bhattacharyya}
#' is \code{FALSE}.
#' @param shrink By default, (if \code{TRUE}), we shrink each covariance
#' matrix with the MDEB covariance matrix estimator. Otherwise, no shrinkage is
#' applied.  Ignored if \code{bhattacharyya} is \code{FALSE}.
#' @param ... Additional arguments passed to the classifier. See details.
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
simdiag.default <- function(x, y, classifier = c("linear", "quadratic"),
                            q = NULL, bhattacharyya = TRUE, fast_svd = TRUE,
                            pool_cov = FALSE, shrink = FALSE, pct = 0.9,
                            tol = 1e-6, ...) {
  x <- as.matrix(x)
  y <- as.factor(y)
  p <- ncol(x)

  classifier <- match.arg(classifier)

  # Creates an object that contains the results returned.
  obj <- list()
  
  # Selects the value for 'q' in one of three ways:
  # 1. If a value for 'q' is provided, we check that it is between 1 and p,
  # inclusively. If not, an error is thrown.
  # 2. Use the Bhattacharyya distance.
  # 3. Use the rank of the sample covariance matrices
  if (bhattacharyya) {
    obj$bhattacharyya <- bhatta_simdiag(x = x, y = y, q = q, pct = pct,
                                        pool_cov = pool_cov, shrink = shrink,
                                        tol = tol)
    obj$q <- obj$bhattacharyya$q
    # Reduce the rank of the simultaneous diagonalizer and the transformed data
    # set.
    simdiag_cov_out <- simdiag_cov(x = x, y = y, q = NULL, fast_svd = fast_svd,
                                 shrink = shrink, tol = tol)

    obj$Q <- simdiag_cov_out$Q[obj$bhattacharyya$dist_rank, ]
    obj$x <- simdiag_cov_out$x[, obj$bhattacharyya$dist_rank]
  } else {
    simdiag_cov_out <- simdiag_cov(x = x, y = y, q = q, fast_svd = fast_svd,
                                   shrink = shrink, tol = tol)

    if (!is.null(q)) {
      q <- as.integer(q)
      if (q < 1 || q > p) {
        stop("The value for 'q' must be between 1 and p, inclusively.")
      }
      obj$q <- q
    } else {
      obj$q <- simdiag_cov_out$q
    }
    
    obj$Q <- simdiag_cov_out$Q[seq_len(obj$q), ]
    obj$x <- simdiag_cov_out$x[, seq_len(obj$q)]
  }

  # Constructs the classifier from the transformed data set using DLDA or DQDA
  # based on the user's specification.
  if (classifier == "linear") {
    obj$classifier <- dlda(x = obj$x, y = y, shrink = shrink, ...)
  } else {
    obj$classifier <- dqda(x = obj$x, y = y, shrink = shrink, ...)
  }

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

#' SimDiag prediction of the class membership of a matrix of new observations.
#'
#' First, we apply the simultaneous diagonalizer constructed from a training
#' data set. Then, we apply the diagonal linear or quadratic classifier based on
#' the user specification to the \code{\link{simdiag}} function.
#' 
#' @rdname simdiag
#' @method predict simdiag
#' @S3method predict simdiag
#' @export
#'
#' @param object object of class \code{simdiag}
#' @param newdata matrix of observations to predict. Each row corresponds to a new observation.
#'
#' @return list predicted class memberships of each row in newdata
predict.simdiag <- function(object, newdata) {
	if (!inherits(object, "simdiag"))  {
		stop("object not of class 'simdiag'")
	}
	if (is.vector(newdata)) {
    newdata <- matrix(newdata, nrow = 1)
  }

  # Transform the test data set with the simultaneous diagonalizer Q.
  newdata <- simdiag_transform(object = object, x = newdata)

  # Predict the class membership of the transformed test data set.
	predict(object = object$classifier, newdata = newdata)
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

#' Plots the Bhattacharyya or Mahalanobis distances of the transformed features
#' between the two classes after simultaneously diagonalizing the sample
#' covariance matrices of each class.
#'
#' The object must be an object of class \code{simdiag}.
#'
#' @keywords internal
#' @param x object of class \code{simdiag} to plot
#' @param ... unused
#' @rdname simdiag
#' @method plot simdiag
#' @S3method plot simdiag
#' @export
plot.simdiag <- function(x, ...) {
  warning("Not yet implemented")
}

