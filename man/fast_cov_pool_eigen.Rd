\name{fast_cov_pool_eigen}
\alias{fast_cov_pool_eigen}
\title{Computes the eigenvalue decomposition of the pooled sample covariance matrix
for the data matrix \code{x} using the Fast Singular Value Decomposition (SVD).}
\usage{
  fast_cov_pool_eigen(x, y, tol = 1e-06)
}
\arguments{
  \item{x}{data matrix with \code{n} observations and
  \code{p} feature vectors}

  \item{y}{class labels for observations (rows) in
  \code{x}}

  \item{tol}{tolerance value below which the singular
  values of \code{x} are considered zero.}
}
\value{
  a list containing: the eigenvectors of the covariance
  matrix MLE the eigenvalues of the covariance matrix MLE
}
\description{
  First, we use the \code{fast.svd} function from the
  \code{corpcor} package to quickly compute the SVD of the
  matrix \code{x}. Then, we calculate the eigenvalues and
  their corresponding eigenvectors from the SVD.
}
\details{
  This function employs a well-known trick for tall data
  (large \code{n}, small \code{p}) and wide data (large
  \code{p}, small \code{n}) to compute the eigenvalue
  decomposition of the covariance matrix estimator from the
  data matrix, \code{x}.

  TODO: Briefly describe the Fast SVD trick. TODO: Describe
  how we use the SVD to compute the eigenvalue
  decomposition.
}

