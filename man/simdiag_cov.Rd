\name{simdiag_cov}
\alias{simdiag_cov}
\title{Simultaneously diagonalizes two real symmetric matrices.}
\usage{
  simdiag_cov(x, y, fast_svd = FALSE, q = NULL,
    reduced_rank = TRUE, shrink = FALSE, tol = 1e-06)
}
\arguments{
  \item{x}{matrix containing the training data. The rows
  are the sample observations, and the columns are the
  features.}

  \item{y}{vector of class labels for each training
  observation}

  \item{fast_svd}{TODO}

  \item{q}{the reduced dimension of the simultaneous
  diagonalizer.}

  \item{reduced_rank}{Should we select \code{q} based on
  the rank of the sample covariance matrices? Ignored if
  \code{q} is given. See details.}

  \item{shrink}{TODO}

  \item{tol}{a value indicating the magnitude below which
  eigenvalues are considered 0.}
}
\value{
  a list containing: the matrix Q that is the \eqn{p \times
  q} simultaneous diagonalizing matrix the matrix x
  consisting of the simultaneously diagonalized data the
  selected value of \code{q}
}
\description{
  This function simultaneously diagonalizes two real
  symmetric matrices. If the matrices are singular, we
  employ a dimension reduction method that simultaneously
  diagonalizes them to a lower dimension, \code{q}.
}
\details{
  The user can manually provides the dimension, \code{q}.
  This is potentially useful if a scree plot is utilized or
  a manual threshold is desired. Otherwise, we select
  \code{q} automatically to be the number of non-zero
  eigenvalues. To determine what constitutes 'non-zero', we
  select all of the eigenvalues larger than the specified
  tolerance, \code{tol}. This is similar to a typical
  automatic approach to Principal Components Analysis.

  The returned simultaneous diagonalizing matrix, \code{Q},
  is constructed such that

  \deqn{Q A Q' = D,} and \deqn{Q B Q' = I_q},

  where \eqn{I_q} denotes the \eqn{q}-dimensional identity
  matrix and \eqn{D} is the diagonal matrix with diagonal
  entries as the characteristic roots (generalized
  eigenvalues) of det(B - \\lambda A) = 0.

  TODO: Discuss briefly that we set 'q' to the number of
  nonzero eigenvalues if the Fast SVD is used.
}

