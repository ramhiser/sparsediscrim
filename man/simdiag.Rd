\name{simdiag}
\alias{print.simdiag}
\alias{simdiag}
\alias{simdiag.default}
\alias{simdiag.formula}
\title{Simultaneously diagonalizes the sample covariance matrices for a data set
generated from \code{K} classes.}
\usage{
  simdiag(x, ...)

  \method{simdiag}{default} (x, y, q = NULL,
    bhattacharyya = TRUE, fast_svd = TRUE, tol = 1e-05)

  \method{simdiag}{formula} (formula, data, q = NULL,
    bhattacharyya = TRUE, fast_svd = TRUE, tol = 1e-05)

  \method{print}{simdiag} (x, ...)
}
\arguments{
  \item{x}{matrix containing the training data. The rows
  are the sample observations, and the columns are the
  features.}

  \item{y}{vector of class labels for each training
  observation}

  \item{q}{the reduced dimension of the simultaneous
  diagonalizer. If \code{NULL} (default), we choose
  \code{q} via the Bhattacharyya distance method, discussed
  below.}

  \item{fast_svd}{logical value that indicates if we should
  utilize the Fast SVD compute the eigenvalue decomposition
  of the sample covariance matrix having larger rank. By
  default, we use this Fast SVD method, which is
  computationally far superior for 'wide data' (large
  \code{p}, small {n}).}

  \item{bhattacharyya}{logical value that indicates if we
  should select the reduced dimension, \code{q}, via our
  proposed Bhattacharyya method. If \code{TRUE}, we choose
  the value for \code{q} to be the dimension, where the
  Bhattacharyya distance between the two classes hardly
  changes for \code{q + 1}. Our technique is similar, in
  spirit, to the 'elbow criterion' in a Principal
  Components Analysis scree plot. The \code{bhattacharyya}
  argument is ignored if a value for \code{q} is
  specified.}

  \item{tol}{a value indicating the magnitude below which
  eigenvalues are considered 0.}

  \item{formula}{formula of the form \code{groups ~ x1 + x2
  + ...} That is, the response is the grouping factor and
  the right hand side specifies the (non-factor)
  discriminators.}

  \item{data}{data frame from which variables specified in
  \code{formula} are preferentially to be taken.}

  \item{x}{object to print}

  \item{...}{unused}
}
\value{
  a list containing: \itemize{ \item \code{Q}: simultaneous
  diagonalizing matrix of size \eqn{p \times q} \item
  \code{q}: the reduced dimension of the data \item
  \code{x}: the transformed data matrix, \code{x}. }
}
\description{
  This function determines the matrix, \code{Q}, that
  simultaneously diagonalizes the sample covariance
  matrices for a data set generated from \code{K} classes.
  If the matrices are singular, we employ a dimension
  reduction method that simultaneously diagonalizes them to
  a lower dimension, \code{q}.

  Summarizes the simultaneous diagonalization information
  for a set of classes.
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

  \deqn{Q \widehat{\\Sigma_k} Q' = D_k}

  for all \eqn{k = 1, \ldots, K}, where
  \eqn{\widehat{\\Sigma_k}} is the maximum likelihood
  estimator (under normality) for the \eqn{k}th class
  covariance matrix, \eqn{\\Sigma_k}, where \eqn{D_k} is a
  diagonal matrix.

  The returned simultaneously-diagonalizing matrix,
  \code{Q}, is of size \eqn{p \times q}.

  Currently, we only allow \eqn{K = 2}. The \eqn{K = 1}
  case is equivalent to PCA. The \eqn{K \ge 2} case is of
  ongoing research; note that a set of \code{K} positive
  semidefinite matrices are simultaneously diagonalizable
  if they are pairwise-commute (i.e. \eqn{A B = B A} for
  any two matrices \eqn{A} and \eqn{B} in the set).
  However, pairwise commuting is a strict assumption that
  we wish to relax.
}
\keyword{internal}

