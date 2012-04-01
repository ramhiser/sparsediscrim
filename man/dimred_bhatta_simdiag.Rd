\name{dimred_bhatta_simdiag}
\alias{dimred_bhatta_simdiag}
\title{Determines the reduced dimension of the simultaneously diagonalized data
based on the Bhattacharyya distance.}
\usage{
  dimred_bhatta_simdiag(x, y, shrink = TRUE, tol = 0.001)
}
\arguments{
  \item{x}{data matrix with \code{n} observations and
  \code{p} feature vectors}

  \item{y}{class labels for observations (rows) in
  \code{x}}

  \item{tol}{tolerance value. If the percentage change in
  the Bhattacharyya distance for a given value of \code{q}
  and \code{q + 1} is less than \code{tol}, then we select
  the value of \code{q} as the optimal value. If there are
  no cases where the percentage change is less than
  \code{tol}, then we select the reduced dimension to be
  the \code{ncol(x)}.}

  \item{shrink}{By default, (if \code{TRUE}), we shrink
  each covariance matrix with the MDEB covariance matrix
  estimator. Otherwise, no shrinkage is applied.}
}
\value{
  a list containing: the selected value of \code{q} the
  estimated Bhattacharyya distance between the two classes
  for each candidate value of \code{q} from 2 to
  \code{ncol(x)}.  the percentage change between the
  Bhattacharyya distances for the sequence of candidate
  values of \code{q}.
}
\description{
  For a data matrix, \code{x}, that has been simultaneously
  diagonalized, we calculate the Bhattacharyya distance
  (divergence) between the two classes given in \code{y}
  for the dimensions ranging from 2 to \code{ncol(x)}. We
  assume that the matrix \code{x} has been simultaneously
  diagonalized.
}
\details{
  We wish to determine the value of \code{q}, such that the
  percentage change in the Bhattacharyya distance for the
  first \code{q} columns and the Bhattacharyya distance for
  the first \code{q + 1} columns is less than \code{tol}.
  The selected value of \code{q} is the dimension to which
  the data matrix, \code{x}, is reduced. If are no cases
  where the percentage change is less than \code{tol}, then
  we select the reduced dimension to be the \code{ncol(x)}.
  Notice that our procedure is similar to that of the
  'elbow criterion' often employed with a scree plot in a
  Principal Components Analysis (PCA). Our approach allows
  us to generalize the PCA dimension reduction technique,
  accordingly.

  By default, we shrink the (diagonal) maximum likelihood
  estimators (MLEs) for the covariance matrices with the
  Minimum Distance Empirical Bayes estimator (MDEB) from
  Srivastava and Kubokawa (2007), which effectively shrinks
  the MLEs towards an identity matrix that is scaled by the
  average of the nonzero eigenvalues.

  TODO: Provide a reference for Fukunaga's text for
  Bhattacharyya distance. TODO: Provide a reference for
  Srivastava and Kubokawa (2007).
}

