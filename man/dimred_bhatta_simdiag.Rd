\name{dimred_bhatta_simdiag}
\alias{dimred_bhatta_simdiag}
\title{Determines the reduced dimension of the simultaneously diagonalized data
based on the Bhattacharyya distance.}
\usage{
  dimred_bhatta_simdiag(x, y, bhatta_prop = 0.9,
    shrink = TRUE)
}
\arguments{
  \item{x}{data matrix with \code{n} observations and
  \code{p} feature vectors}

  \item{y}{class labels for observations (rows) in
  \code{x}}

  \item{bhatta_prop}{threshold value for the maximum
  cumulative proportion of the Bhattacharyya distance. We
  use the threshold to determine the reduced dimension
  \code{q}. When paired with the SimDiag method, the
  Bhattacharyya approach to dimension reduction generalizes
  the scree plot idea that is often utilized with Principal
  Components Analysis.}

  \item{shrink}{By default, (if \code{TRUE}), we shrink
  each covariance matrix with the MDEB covariance matrix
  estimator. Otherwise, no shrinkage is applied.}
}
\value{
  a list containing: \itemize{ \item \code{q}: the reduced
  dimension determined via \code{bhatta_pct} \item
  \code{bhatta_dist}: the estimated Bhattacharyya distance
  between the two classes for each dimension in the
  transformed space.  \item \code{bhatta_dist_rank}: the
  indices of the columns of \code{x} corresponding to the
  \code{q} largest Bhattacharyya distances.  \item
  \code{bhatta_cumprop}: the cumulative proportion of the
  sorted Bhattacharyya distances for each column given in
  \code{x}. }
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
  We wish to determine the \code{q} largest Bhattacharyya
  distances from the \code{p} features in \code{x}. We
  assume that the class covariance matrices for both
  classes given in \code{y} are diagonal. The selected
  value of \code{q} is the dimension to which the data
  matrix, \code{x}, is reduced. Notice that our procedure
  is similar to the reduced dimension selection technique
  often employed with a scree plot in a Principal
  Components Analysis (PCA). Our approach allows us to
  generalize the PCA dimension reduction technique,
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

