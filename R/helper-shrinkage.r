#' Applies covariance matrix shrinkage to each class for the DLDA and DQDA
#' classifiers.
#'
#' For each class covariance matrix (which is assumed to be diagonal),
#' we shrink the diagonal matrix towards a diagonal matrix, \eqn{D_k}.
#' We assume that each sample covariance matrix, \eqn{\widehat{\Sigma}_k}, is
#' diagonal for \eqn{k = 1, \ldots, K}, where \eqn{K} is the number of classes.
#' We wish to shrink each covariance matrix towards a diagonal matrix,
#' \eqn{D_k}. The matrix \eqn{D_k} is given as the \eqn{k}th list element in the
#' list, given as an argument, \code{shrinkage}.
#'
#' @export
#' @param obj object of type \code{dlda} or \code{dqda}
#' @param shrinkage list of length \eqn{K}. The \eqn{k}th list element is a
#' diagonal matrix towards which the \eqn{k}th sample covariance matrix is
#' shrunken.
#' @return list of length \code{K} with MDEB shrinkage matrices
diag_shrinkage <- function(obj, diag_mat = rep(1, p), pool = FALSE) {
  if (pool) {
    obj$var_pool <- obj$var_pool + shrink
    obj$shrinkage <- shrink
  } else {
    obj$est <- mapply(function(class_est, shrink) {
      class_est$var <- class_est$var + shrink
      class_est$shrinkage <- shrink
      class_est
    }, obj$est, shrinkage)
  }
  obj
}

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
#' @param diag_mat vector of length \code{p} that is the diagonal matrix towards
#' which each class' covariance matrix is shrunken
#' @return list of length \code{K} with MDEB shrinkage matrices
mdeb_shrinkage <- function(obj, diag_mat = rep(1, p)) {
  p <- obj$p
  lapply(obj$est, function(class_est) {
    shrink_coeff <- sum(class_est$var)
    shrink_coeff * diag_mat
  })
}
