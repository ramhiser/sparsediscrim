\name{simdiag_fastsvd}
\alias{simdiag_fastsvd}
\title{Simultaneously diagonalizes the sample covariance matrices for a data set
generated from 2 classes with the Fast SVD method.}
\usage{
  simdiag_fastsvd(x, y, calc_ranks = FALSE, q = NULL,
    tol = 1e-06)
}
\arguments{
  \item{x}{matrix containing the training data. The rows
  are the sample observations, and the columns are the
  features.}

  \item{y}{vector of class labels for each training
  observation}

  \item{calc_ranks}{logical value. We consider the first
  class to be the class with sample covariance of larger
  rank. By default, if \code{calc_ranks} is \code{FALSE},
  we consider the first class to be the one with a larger
  sample size.in terms of the larger. Otherwise, if
  \code{TRUE}, we calculate the ranks numerically as the
  number of nonzero eigenvalues and select the first class
  to have larger rank. Note that this option adds
  additional computations. If multicollinearity is
  suspected and if the sample sizes are approximately
  equal, we recommend that \code{calc_ranks} be set to
  \code{TRUE}.}

  \item{q}{the reduced dimension of the simultaneous
  diagonalizer. By default, \code{q} is selected to be the
  rank of the sample covariance matrix for the first
  class.}

  \item{tol}{a value indicating the magnitude below which
  eigenvalues are considered 0.}
}
\value{
  matrix of size \eqn{q \times p} that simultaneously
  diagonalizes the sample covariance matrices (MLEs) for
  each class
}
\description{
  This function determines the matrix, \code{Q}, that
  simultaneously diagonalizes the sample covariance
  matrices for a data set generated from 2 classes. If the
  matrices are singular, we employ a dimension reduction
  method that simultaneously diagonalizes them to a lower
  dimension, \code{q}. This function simultaneously
  diagonalizes two real symmetric matrices. If the matrices
  are singular, we employ a dimension reduction method that
  simultaneously diagonalizes them to a lower dimension,
  \code{q}.
}
\details{
  The user can manually provides the dimension, \code{q}.
  This is potentially useful if a scree plot is utilized or
  a manual threshold is desired. Otherwise, we select
  \code{q} automatically to be the number of non-zero
  eigenvalues. To determine what constitutes 'non-zero', we
  select all of the eigenvalues larger than the specified
  tolerance, \code{tol}. This is similar to a typical
  automatic approach to Principal Components Analysis
  (PCA).

  The returned simultaneous diagonalizing matrix, \code{Q},
  is constructed such that

  \deqn{Q \\Sigma_1 Q' = I_q,} and \deqn{Q \\Sigma_2 Q' =
  D},

  where \eqn{I_q} denotes the \eqn{q}-dimensional identity
  matrix and \eqn{D} is the diagonal matrix with diagonal
  entries as the characteristic roots (generalized
  eigenvalues) of det(\\Sigma_1 - \\lambda \\Sigma_2) = 0.

  We use the 'Fast Singular Value Decomposition (SVD)'
  approach to determine the eigenvectors and eigenvalues of
  the class with the larger rank. This technique provides
  significant improvements in computational runtime when
  the number of observations number is much less than the
  number of features. That is, we gain significantly in
  terms of runtime if the number of rows of \code{x} is
  much less than the number of columns of \code{x}.

  We remark that the rank estimation (i.e. choosing the
  value of \code{q} when it is not manually specified) can
  be subject to numerical instabilities if the value for
  \code{tol} is too small. In fact, the ranks of the
  covariance matrices can be overestimated, and therefore,
  the value for \code{q} can be too larger; this issue
  yields a linear transformation, \code{Q}, that projects
  the data to a dimension larger than the actual rank. The
  result is that sample covariance matrices might not be
  fully simultaneously diagonalized. Our approach to
  overcome the issue is to set \code{q} to \code{max(n1,
  n2) - 1} if \code{q >= max(n1, n2)} because the ranks of
  the class sample covariance matrices are bounded above by
  the class sample size minus 1.
}

