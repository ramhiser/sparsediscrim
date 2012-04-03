#' Computes the MDEB shrinkage matrices for each of \code{K} classes.
#'
#' For each class covariance matrix (which is assumed to be diagonal),
#' we compute the diagonal matrix towards each class covariance matrix is
#' shrunken. The resulting matrix is scaled by the MDEB estimator for the
#' shrinkage parameter, \eqn{\gamma}.
#'
#' We assume that each sample covariance matrix, \eqn{\widehat{\Sigma}_k}, is
#' diagonal for \eqn{k = 1, \ldots, K}, where \eqn{K} is the number of classes.
#' We wish to shrink each covariance matrix towards a common diagonal matrix,
#' \eqn{D}. The matrix \eqn{D} is given as the argument \code{diag_mat}.
#'
#' Ultimately, we construct the following shrunken covariance matrix estimator:
#' \deqn{\widehat{\Sigma}_k(\gamma) = \widehat{\Sigma}_k + \gamma_k D}.
#' The MDEB estimator for \eqn{\gamma_k} is given by:
#' \deqn{\gamma_k = \text{trace}(\widehat{\Sigma}_k) / \min{n_k, p}},
#' which is approximately equal to the average of the nonzero eigenvalues.
#'
#' Explicitly, this function returns a list, where each list element contains
#' \eqn{\gamma_k D}, where \eqn{\gamma_k} is given by the above MDEB estimator.
#'
#' TODO: Cite the Srivastava and Kubokawa paper.
#' 
#' @export
#' @param obj object of type \code{dlda} or \code{dqda}
#' @param pool logical value. If \code{TRUE}, the population covariance matrices
#' are considered equal, and the MDEB shrinkage is applied to the (diagonal)
#' pooled sample covariance matrix in the \code{var_pool} vector. If \code{FALSE}
#' (default), the MDEB shrinkage estimator is computed for each class' (diagonal)
#' covariance matrix in its corresponding vector, \code{var}.
#' @return list of length \code{K} with MDEB shrinkage matrices
mdeb_shrinkage <- function(obj, pool = FALSE) {
  if (pool) {
    obj$shrinkage <- with(obj, sum(var_pool) + min(p, N))
    obj$var_pool <- with(obj, var_pool + shrinkage)
    
  } else {
    obj$est <- lapply(obj$est, function(class_est) {
      class_est$shrinkage <- with(class_est, sum(var) / min(obj$p, n))
      class_est$var <- with(class_est, var + shrinkage)
      class_est
    })
  }

  obj
}
