\name{cov_block_autocorrelation}
\alias{cov_block_autocorrelation}
\title{Generates a \eqn{p \times p} block-diagonal covariance matrix with
autocorrelated blocks.}
\usage{
  cov_block_autocorrelation(num_blocks, block_size, rho,
    sigma2 = 1)
}
\arguments{
  \item{num_blocks}{the number of blocks in the covariance
  matrix}

  \item{block_size}{the size of each square block within
  the covariance matrix}

  \item{rho}{the autocorrelation parameter. Must be less
  than 1 in absolute value.}

  \item{sigma2}{the variance of each feature}
}
\value{
  autocorrelated covariance matrix
}
\description{
  This function generates a \eqn{p \times p} covariance
  matrix with autocorrelated blocks. The autocorrelation
  parameter is \code{rho}. There are \code{num_blocks}
  blocks each with size, \code{block_size}. The variance,
  \code{sigma2}, is constant for each feature and defaulted
  to 1.
}
\details{
  The autocorrelated covariance matrix is defined as: TODO:
  Give definition.

  The value of \code{rho} must be such that \eqn{\abs{\rho}
  < 1} to ensure that the covariance matrix is positive
  definite.

  The size of the resulting matrix is \eqn{p \times p},
  where \code{p = num_blocks * block_size}.
}

