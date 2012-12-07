\name{simdiag_cov}
\alias{simdiag_cov}
\title{Simultaneously diagonalizes two real symmetric matrices.}
\usage{
  simdiag_cov(A, B, q = NULL, tol = 1e-06)
}
\arguments{
  \item{A}{p x p real symmetric matrix}

  \item{B}{p x p real symmetric matrix}

  \item{q}{the reduced dimension of the simultaneous
  diagonalizer.}

  \item{tol}{a value indicating the magnitude below which
  eigenvalues are considered 0.}
}
\value{
  a list containing: the matrix Q that is the \eqn{p \times
  q} simultaneous diagonalizing matrix the selected value
  of \code{q}
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
}

