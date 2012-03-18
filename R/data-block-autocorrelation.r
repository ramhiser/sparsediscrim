#' Generates data from \code{K} multivariate normal data populations, where each
#' population (class) has a covariance matrix consisting of block-diagonal
#' autocorrelation matrices.
#'
#' This function generates \code{K} multivariate normal data sets, where each
#' data set is generated with a constant mean vector and a covariance matrix
#' consisting of block-diagonal autocorrelation matrices. The data are returned
#' as a single matrix, \code{x}, along with a vector of class labels, \code{y},
#' that indicates class membership.
#' 
#' For simplicity, we assume that a class mean vector is constant for each
#' feature. That is, we assume that the mean vector of the \eqn{k}th class is
#' \eqn{c_k * j_p}, where \eqn{j_p} is a \eqn{p \times 1} vector of ones and
#' \eqn{c_k} is a real scalar.
#' 
#' TODO: Define the covariance matrix, \eqn{\Sigma_k}.
#'
#' The number of classes, \code{K}, is determined with lazy evaluation as the
#' \code{length} of \code{n}.
#'
#' @export
#' @param n vector of the sample sizes of each class. The \code{length} of
#' \code{n} determines the number of classes, \code{K}.
#' @param block_size TODO
#' @param num_blocks TODO
#' @param rho vector of the values of the autocorrelation parameter for each
#' class covariance matrix. Must equal the \code{length} of \code{n} (i.e. equal
#' to \code{K}).
#' @param sigma2 vector of the variance coefficients for each class covariance
#' matrix. Must equal the \code{length} of \code{n} (i.e. equal to \code{K}).
#' @return named list with elements:
#' \itemize{
#'   \item \code{x}: matrix of observations with \code{sum(n)} rows and \code{p}
#' columns
#'   \item \code{y}: vector of class labels that indicates class membership for
#' each observation (row) in \code{x}.
#' }
#' @examples
#' # Generates data from K = 3 classes.
#' data <- generate_blockdiag(n = 3:5, p = 5, rho = seq(.1, .9, length = 3),
#'                             mu = c(0, 3, -2))
#' data$x
#' data$y
#' 
#' # Generates data from K = 4 classes. Notice that we use specify a variance.
#' data <- generate_blockdiag(n = 3:6, p = 4, rho = seq(0, .9, length = 4),
#'                             mu = c(0, 3, -2, 6), sigma2 = 1:4)
#' data$x
#' data$y
generate_blockdiag <- function(n, block_size, num_blocks, rho, mu,
                               sigma2 = rep(1, K)) {
  require('mvtnorm')
  require('bdsmatrix')
  n <- as.integer(n)
  block_size <- as.integer(block_size)
  num_blocks <- as.integer(num_blocks)
  p <- block_size * num_blocks
  rho <- as.numeric(rho)
  mu <- as.numeric(mu)

  K <- length(n)

  if (length(rho) != K) {
    stop("The length of 'rho' must equal the length of 'n'.")
  } else if(length(mu) != K) {
    stop("The length of 'mu' must equal the length of 'n'.")
  }
  x <- lapply(seq_len(K), function(k) {
    mvtnorm:::rmvnorm(n = n[k], mean = rep(mu[k], p),
                      sigma = cov_block_autocorrelation(num_blocks = num_blocks,
                        block_size = block_size, rho = rho[k],
                        sigma2 = sigma2[k]))
  })
  x <- do.call(rbind, x)
  y <- factor(rep(seq_along(n), n))
  
  list(x = x, y = y)
}

#' Generates a \eqn{p \times p} block-diagonal covariance matrix with
#' autocorrelated blocks.
#'
#' This function generates a \eqn{p \times p} covariance matrix with
#' autocorrelated blocks. The autocorrelation parameter is \code{rho}.
#' There are \code{num_blocks} blocks each with size, \code{block_size}.
#' The variance, \code{sigma2}, is constant for each feature and defaulted to 1.
#'
#' The autocorrelated covariance matrix is defined as:
#' TODO: Give definition.
#'
#' The value of \code{rho} must be such that \eqn{\abs{\rho} < 1} to ensure that
#' the covariance matrix is positive definite.
#'
#' The size of the resulting matrix is \eqn{p \times p}, where
#' \code{p = num_blocks * block_size}.
#'
#' @param num_blocks the number of blocks in the covariance matrix
#' @param block_size the size of each square block within the covariance matrix
#' @param rho the autocorrelation parameter. Must be less than 1 in absolute
#' value.
#' @param sigma2 the variance of each feature
#' @return autocorrelated covariance matrix
cov_block_autocorrelation <- function(num_blocks, block_size, rho, sigma2 = 1) {
  require('bdsmatrix')
  cov_block <- diagdiscrim:::cov_autocorrelation(p = block_size, rho = rho,
                                                 sigma2 = sigma2)
  cov_block <- as.vector(cov_block)
  as.matrix(bdsmatrix:::bdsmatrix(blocksize = rep(block_size, num_blocks),
                        blocks = replicate(num_blocks, cov_block)))
}

#' Generates a \eqn{p \times p} autocorrelated covariance matrix
#'
#' This function generates a \eqn{p \times p} autocorrelated covariance matrix
#' with autocorrelation parameter, \code{rho}. The variance, \code{sigma2}, is constant for each
#' feature and defaulted to 1.
#'
#' The autocorrelated covariance matrix is defined as:
#' TODO: Give definition.
#'
#' The value of \code{rho} must be such that \eqn{\abs{\rho} < 1} to ensure that
#' the covariance matrix is positive definite.
#'
#' @param p the size of the covariance matrix
#' @param rho the autocorrelation parameter. Must be less than 1 in absolute
#' value.
#' @param sigma2 the variance of each feature
#' @return autocorrelated covariance matrix
cov_autocorrelation <- function(p, rho, sigma2 = 1) {
  p <- as.integer(p)
  rho <- as.numeric(rho)
  sigma2 <- as.numeric(sigma2)
  
  if (abs(rho) >= 1) {
    stop("The value of 'rho' must be between (1 - p)^(-1) and 1, exclusively.")
  }
  if (sigma2 <= 0) {
    stop("The value of 'sigma2' must be positive.")
  }
  x <- diag(p)
  sigma2 * rho^abs(row(x) - col(x))
}
