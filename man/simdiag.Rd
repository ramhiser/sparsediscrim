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
    fast_svd = TRUE, bhattacharyya = TRUE, tol = 1e-05)

  \method{simdiag}{formula} (formula, data, q = NULL,
    method = c("rank", "bhattacharyya"), tol = 1e-06)

  \method{print}{simdiag} (x, ...)
}
\arguments{
  \item{x}{matrix containing the training data. The rows
  are the sample observations, and the columns are the
  features.}

  \item{y}{vector of class labels for each training
  observation}

  \item{q}{the reduced dimension of the simultaneous
  diagonalizer.}

  \item{method}{the method to automatically select the
  value of \code{q}. Ignored if \code{q} is given.}

  \item{simdiag_data}{Should the data matrix, \code{x}, be
  transformed and returned? If TRUE (default), we return an
  updated matrix, \code{x}. For large data sets, this
  option can yield a very large \code{simdiag} object and
  should possibly be calculated manually.}

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
  \code{TRUE}. Only used if the \code{method} is selected
  to be \code{fast_svd}; this is the default option.}

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
  \code{x}: the transformed data matrix, \code{x}.
  Returned, if \code{simdiag_data} is TRUE. Otherwise,
  \code{x} is not returned. }
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

