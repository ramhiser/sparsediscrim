\name{bhattacharyya}
\alias{bhattacharyya}
\title{Estimates the Bhattacharyya Distance between Two Multivariate Normal
Populations}
\usage{
  bhattacharyya(x, y, diag = FALSE, pool_cov = FALSE,
    shrink = FALSE)
}
\arguments{
  \item{x}{data matrix with \code{n} observations and
  \code{p} feature vectors}

  \item{y}{class labels for observations (rows) in
  \code{x}}

  \item{diag}{logical. Should we assume that the covariance
  matrices are diagonal? If so, this greatly simplifies and
  expedites the distance calculation}
}
\value{
  the estimated Bhattacharyya distance between the two
  classes given in \code{y}.
}
\description{
  For a data matrix, \code{x}, we calculate the
  Bhattacharyya distance (divergence) between the two
  classes given in \code{y}.
}
\details{
  An excellent overview of the Bhattacharyya distance can
  be found in Fukunaga (1990).
}
\references{
  Fukunaga, Keinosuke (1990). Introduction to Statistical
  Pattern Recognition. Academic Press Inc., 2nd edition,
  1990. \url{http://amzn.to/Ke9m6l}.
}

