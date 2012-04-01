\name{bhatta_simdiag}
\alias{bhatta_simdiag}
\title{Bhattacharyya Distance for Simultaneously Diagonalized Data}
\usage{
  bhatta_simdiag(x, y, q = NULL, shrink = TRUE)
}
\arguments{
  \item{x}{data matrix with \code{n} observations and
  \code{p} feature vectors}

  \item{y}{class labels for observations (rows) in
  \code{x}}

  \item{q}{the number of columns of \code{x} to use in
  calculating the Bhattacharyya distance. By default, we
  use \code{ncol(x)}.}

  \item{shrink}{By default, (if \code{TRUE}), we shrink
  each covariance matrix with the MDEB covariance matrix
  estimator. Otherwise, no shrinkage is applied.}
}
\value{
  the estimated Bhattacharyya distance between the two
  classes given in \code{y}.
}
\description{
  For a data matrix, \code{x}, that has been simultaneously
  diagonalized, we calculate the Bhattacharyya distance
  (divergence) between the two classes given in \code{y}.
}
\details{
  We assume that the matrix \code{x} has been
  simultaneously diagonalized. If simultaneous
  diagonalization has not been employed, this function, by
  default, will calculate the Bhattacharyya distance
  between the two classes under the assumption that each
  class is multivariate normal with diagonal covariance
  matrices. However, in this case any choice for \code{q},
  other than the default, will not be applicable.

  By default, we shrink the (diagonal) maximum likelihood
  estimators (MLEs) for the covariance matrices with the
  Minimum Distance Empirical Bayes estimator (MDEB) from
  Srivastava and Kubokawa (2007), which effectively shrinks
  the MLEs towards an identity matrix that is scaled by the
  average of the nonzero eigenvalues.

  The value of \code{q} determines the number of columns
  that will be used to calculate the Bhattacharyya
  distance. This approach is useful in determining an
  'optimal' value for \code{q} over a range of candidate
  values.

  TODO: Provide a reference for Fukunaga's text for
  Bhattacharyya distance. TODO: Provide a reference for
  Srivastava and Kubokawa (2007).
}

