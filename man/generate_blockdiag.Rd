\name{generate_blockdiag}
\alias{generate_blockdiag}
\title{Generates data from \code{K} multivariate normal data populations, where each
population (class) has a covariance matrix consisting of block-diagonal
autocorrelation matrices.}
\usage{
  generate_blockdiag(n, block_size, num_blocks, rho, mu,
    sigma2 = rep(1, K))
}
\arguments{
  \item{n}{vector of the sample sizes of each class. The
  \code{length} of \code{n} determines the number of
  classes, \code{K}.}

  \item{block_size}{TODO}

  \item{num_blocks}{TODO}

  \item{rho}{vector of the values of the autocorrelation
  parameter for each class covariance matrix. Must equal
  the \code{length} of \code{n} (i.e. equal to \code{K}).}

  \item{sigma2}{vector of the variance coefficients for
  each class covariance matrix. Must equal the
  \code{length} of \code{n} (i.e. equal to \code{K}).}
}
\value{
  named list with elements: \itemize{ \item \code{x}:
  matrix of observations with \code{sum(n)} rows and
  \code{p} columns \item \code{y}: vector of class labels
  that indicates class membership for each observation
  (row) in \code{x}. }
}
\description{
  This function generates \code{K} multivariate normal data
  sets, where each data set is generated with a constant
  mean vector and a covariance matrix consisting of
  block-diagonal autocorrelation matrices. The data are
  returned as a single matrix, \code{x}, along with a
  vector of class labels, \code{y}, that indicates class
  membership.
}
\details{
  For simplicity, we assume that a class mean vector is
  constant for each feature. That is, we assume that the
  mean vector of the \eqn{k}th class is \eqn{c_k * j_p},
  where \eqn{j_p} is a \eqn{p \times 1} vector of ones and
  \eqn{c_k} is a real scalar.

  TODO: Define the covariance matrix, \eqn{\Sigma_k}.

  The number of classes, \code{K}, is determined with lazy
  evaluation as the \code{length} of \code{n}.
}
\examples{
# Generates data from K = 3 classes.
data <- generate_blockdiag(n = 3:5, p = 5, rho = seq(.1, .9, length = 3),
                            mu = c(0, 3, -2))
data$x
data$y

# Generates data from K = 4 classes. Notice that we use specify a variance.
data <- generate_blockdiag(n = 3:6, p = 4, rho = seq(0, .9, length = 4),
                            mu = c(0, 3, -2, 6), sigma2 = 1:4)
data$x
data$y
}

