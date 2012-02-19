#' Simultaneously diagonalizes two real symmetric matrices.
#'
#' This function simultaneously diagonalizes two real symmetric matrices.
#' If the matrices are singular, we employ a dimension reduction
#' method that simultaneously diagonalizes them in a lower dimension.
#'
#' If the user selects the dimension reduction option, the user may
#' specify a value of k. This is potentially useful if a scree plot
#' is utilized or a manual threshold is desired. Otherwise, we
#' select k automatically to be the number of non-zero eigenvalues.
#' To determine what constitutes 'non-zero', we select all of the eigenvalues
#' larger than the specified tolerance, 'tol'.
#' 
#' If the user does not wish to use the dimension reduction method,
#' then k is automatically computed as the number of eigenvalues of B.
#' Hence, in this case, k = p.
#'
#' The returned simultaneous diagonalizing matrix, Q, is constructed such that
#'
#' Q A Q' = D, and
#' Q B Q' = I_k,
#'
#' where I_k denotes the k-dimensional identity matrix and Lambda is the 
#' diagonal matrix with diagonal entries as the characteristic roots
#' (generalized eigenvalues) of det(B - lambda A) = 0.
#'
#' @export
#' @param A p x p real symmetric matrix
#' @param B p x p real symmetric matrix

#' @param q the reduced dimension of the simultaneous diagonalizer.
#' @param method the method to automatically select the value of \code{q}.
#' Ignored if \code{q} is given.
#' @param simdiag_data Should the data matrix, \code{x}, be transformed and
#' returned? If TRUE (default), we return an updated matrix, \code{x}. For large
#' data sets, this option can yield a very large \code{simdiag} object.
#' @param tol a value indicating the magnitude below which eigenvalues are
#' considered 0.
#' @return a list containing:
#' \itemize{
#'   \item \code{Q}: simultaneous diagonalizing matrix of size p x q
#'	 \item \code{q}: the reduced dimension of the data
#'   \item \code{x}: the transformed data matrix, \code{x}. Returned, if
#'   simdiag_data is TRUE. Otherwise, \code{x} is not returned.
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

  # In future work, we intend to allow for more than the K = 2 case. Hence, we
  # have delegated the actual simultaneous diagonalization to the function,
  # 'simdiag2'.
  # TODO: Pass 'method' to simdiag2. For now, we are only using the rank
  # estimation to select 'q'.
  obj <- diagdiscrim:::simdiag2(x = x, y = y, q = q, tol = tol)

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

