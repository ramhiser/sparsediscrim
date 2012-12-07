\name{bhatta_simdiag}
\alias{bhatta_simdiag}
\title{Dimension reduction of two simultaneously diagonalized populations based on
the Bhattacharyya distance.}
\usage{
  bhatta_simdiag(x, y, pct = 0.9, pool_cov = FALSE,
    shrink = FALSE, tol = 1e-06)
}
\arguments{
  \item{x}{data matrix with \code{n} observations and
  \code{p} feature vectors}

  \item{y}{class labels for observations (rows) in
  \code{x}}

  \item{pct}{threshold value for the maximum cumulative
  proportion of the Bhattacharyya distance. We use the
  threshold to determine the reduced dimension \code{q}.
  When paired with the SimDiag method, the Bhattacharyya
  approach to dimension reduction generalizes the scree
  plot idea that is often utilized with Principal
  Components Analysis.}

  \item{shrink}{If \code{TRUE}, we shrink each covariance
  matrix with the MDEB covariance matrix estimator.
  Otherwise, no shrinkage is applied.}

  \item{tol}{a value indicating the magnitude below which
  eigenvalues are considered 0.}
}
\value{
  a list containing: \itemize{ \item \code{q}: the reduced
  dimension determined via \code{bhatta_pct} \item
  \code{dist}: the estimated Bhattacharyya distance between
  the two classes for each dimension in the transformed
  space.  \item \code{dist_rank}: the indices of the
  columns of \code{x} corresponding to the \code{q} largest
  Bhattacharyya distances.  \item \code{cumprop}: the
  cumulative proportion of the sorted Bhattacharyya
  distances for each column given in \code{x}. }
}
\description{
  For a data matrix, \code{x}, we simultaneously
  diagonalize the two classes given in \code{y}. Then, we
  calculate the Bhattacharyya distance (divergence) between
  the two populations for the dimensions ranging from 2 to
  \code{ncol(x)}.
}
\details{
  We wish to determine the \code{q} largest Bhattacharyya
  distances from the \code{p} features in \code{x}. The
  selected value of \code{q} is the dimension to which the
  data matrix, \code{x}, is reduced. Notice that our
  procedure is similar to the reduced dimension selection
  technique often employed with a scree plot in a Principal
  Components Analysis (PCA). Our approach generalizes the
  PCA dimension reduction technique, accordingly.

  An excellent overview of the Bhattacharyya distance can
  be found in Fukunaga (1990).

  We allow shrinkage to be applied to the (diagonal)
  maximum likelihood estimators (MLEs) for the covariance
  matrices with the Minimum Distance Empirical Bayes
  estimator (MDEB) from Srivastava and Kubokawa (2007),
  which effectively shrinks the MLEs towards an identity
  matrix that is scaled by the average of the nonzero
  eigenvalues.
}
\references{
  Fukunaga, Keinosuke (1990). Introduction to Statistical
  Pattern Recognition. Academic Press Inc., 2nd edition,
  1990. \url{http://amzn.to/Ke9m6l}.

  Srivastava, M. S. and Kubokawa, T. (2007). Comparison of
  discrimination methods for high dimensional data. J.
  Japan Statist. Soc., 37(1), 123–134.
}

