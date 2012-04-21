#' Simultaneously diagonalizes the sample covariance matrices for a data set
#' generated from \code{K} classes with one of several algorithmic approaches.
#'
#' This function determines the matrix, \code{Q}, that (nearly) simultaneously
#' diagonalizes the sample covariance matrices for a data set generated from
#' \code{K} classes.
#'
#' The returned simultaneous diagonalizing matrix, \code{Q}, is constructed such
#' that
#'
#' \deqn{Q \widehat{\\Sigma_k} Q' \approx D_k}
#'
#' for all \eqn{k = 1, \ldots, K}, where \eqn{\widehat{\\Sigma_k}} is the maximum
#' likelihood estimator (under normality) for the \eqn{k}th class covariance
#' matrix, \eqn{\\Sigma_k}, where \eqn{D_k} is a diagonal matrix.
#'
#' The returned simultaneously-diagonalizing matrix, \code{Q}, is of size
#' \eqn{p \times p}.
#'
#' TODO: Briefly discuss each algorithmic approach.
#' TODO: Provide reference to implementations in the 'jointDiag' package.
#' @export
#' @param x matrix containing the training data. The rows are the sample
#' observations, and the columns are the features.
#' @param y vector of class labels for each training observation
#' @param method character. The SimDiag algorithm to employ. See Details for more
#' information.
#' @param shrink logical. If \code{TRUE}, each covariance matrix is shrunken
#' toward the a scaled identity matrix to obtain a ridge-like covariance matrix
#' estimator. We use the MDEB estimator to determine the shrinkage term. By
#' default, the argument is \code{FALSE}, and no shrinkage is applied.
#' @param tol a tolerance value below which covergence is reached for the
#' specified algorithm.
#' @return a list containing:
#' \itemize{
#'   \item \code{Q}: simultaneous diagonalizing matrix of size \eqn{p \times p}
#'   \item \code{x}: the transformed data matrix, \code{x}.
#'   \item \code{method}: string with the name of the employed SimDiag algorithm
#' }
#' @examples
#' # We use the well-known Iris data set to demonstrate four algorithmic
#' # approaches to SimDiag.
#' x <- iris[, -5]
#' y <- iris[, 5]
#'
#' # Each of the methods return a list. We are interested in the SimDiag'd
#' # matrix, 'x'.
#' asfari_x <- simdiag_algorithmic(x = x, y = y, method = "asfari")$x
#' jedi_x <- simdiag_algorithmic(x = x, y = y, method = "jedi")$x
#' qdiag_x <- simdiag_algorithmic(x = x, y = y, method = "qdiag")$x
#' ffdiag_x <- simdiag_algorithmic(x = x, y = y, method = "ffdiag")$x
#'
#' # Now, we should have that the covariance matrices for each class are nearly
#' # simultaneously diagonalized. Our goal is to verify that. To do so, we start
#' # off by computing each class's sample covariance matrix MLE for the
#' # SimDiag'd data. We also do the same for the original data.
#' orig_cov <- cov_list(x, y)
#' asfari_cov <- cov_list(asfari_x, y)
#' jedi_cov <- cov_list(jedi_x, y)
#' qdiag_cov <- cov_list(qdiag_x, y)
#' ffdiag_cov <- cov_list(ffdiag_x, y)
#'
#' # To determine how well each SimDiag algorithm performed, we compute the
#' # Frobenius norm between each SimDiag'd covariance matrix and the
#' # corresponding matrix with only its diagonal entries. If the matrix were
#' # fully SimDiag'd, then the compute Frobenius norm should be equal to 0.
#' frob_diag <- function(x) {
#'   norm(x - diag(diag(x)), type = "F")
#' }
#'
#' # We determine that the method that performed best is the one that has the
#' # minimum mean of the Frobenius norms. Note that this is a typical objective
#' # function for SimDiag methods. Note that the \code{qdiag} method achieves
#' # the minimum value.
#' mean(sapply(orig_cov, frob_diag))
#' mean(sapply(asfari_cov, frob_diag))
#' mean(sapply(jedi_cov, frob_diag))
#' mean(sapply(qdiag_cov, frob_diag))
#' mean(sapply(ffdiag_cov, frob_diag))
simdiag_algorithmic <- function(x, y,
                                method = c("asfari", "jade", "jedi", "qdiag",
                                  "ffdiag", "jadiag", "uwedge"), shrink = FALSE,
                                tol = .Machine$double.eps, max_iter = 250) {
  require('jointDiag')
  x <- as.matrix(x)
  y <- as.factor(y)
  method <- match.arg(method)

  # Constructs a list of each class's sample covariance matrix (the MLEs).
  l_covs <- cov_list(x, y, shrink = shrink)

	if(method == "asfari") {
    # The current Asfari implementation requires that the sample covariance
    # matrices be concatenated into a single matrix of size, 'p x Kp', where 'p'
    # is the number of features and 'K' is the number of classes.
		asfari_cov <- do.call(cbind, l_covs)
		Q <- LUJID(asfari_cov, mode = 'B', ERR = tol, RBALANCE = 3, ITER = max_iter)
	} else if(method == "jade") {
		stop("JADE algorithm has not been implemented yet.")
	} else { # The implemented methods from the 'jointDiag' package.
    # The implementations of the algorithmic approaches to SimDiag in the
    # 'jointDiag' package require that the sample covariance matrices be
    # formulated into a 3-dimensional array of size 'p x p x K', where 'p' is the
    # number of features and K is the number of classes. To do this, we use the
    # 'simplify2array'
    jointDiag_cov <- simplify2array(l_covs, higher = TRUE)
    
    # The 'ajd' function acts as the 'jointDiag' package's wrapper function for
    # all of its SimDiag implementations.
		ajd_out <- ajd(jointDiag_cov, method = method, eps = tol, itermax = max_iter)

    # For some crazy reason, all of the methods return the matrix, 'B', as the
    # SimDiag transform matrix except for the Jedi method, which returns the
    # matrix, 'A'.
    if(method == "jedi") {
      Q <- ajd_out$A
    } else {
      Q <- ajd_out$B
    }
	}

	list(Q = Q, x = x %*% Q, method = method)
}

