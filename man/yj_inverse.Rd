\name{yj_inverse}
\alias{yj_inverse}
\title{Computes the inverse of Yeo-Johnson transformed random variate(s)}
\usage{
  yj_inverse(y, lambda, tol = sqrt(.Machine$double.eps))
}
\arguments{
  \item{y}{vector of random variate(s)}

  \item{lambda}{numeric Yeo-Johnson power transformation
  parameter value}

  \item{tol}{numeric tolerance value that determines
  numerical equality}
}
\value{
  the inverse of the Yeo-Johnson random variate, \code{y}
}
\description{
  For a random variate that has been transformed with the
  Yeo-Johnson transformation to near-normality, we invert
  the random variate back to its original scale in this
  function. The inversion requires the transformation value
  \code{lambda} applied initially.
}
\details{
  The YJ transformation to near-normality is a
  generalization of the Box-Cox (BC) transformation that
  allows for observations to be unbounded on the real line,
  whereas the BC transformation requires that observations
  be nonnegative.

  The value for \code{lambda} is typically obtained using
  the \code{\link[car]{powerTransform}} function. The
  transformed random variate is usually obtained from the
  \code{\link[car]{yjPower}} function.
}
\references{
  Yeo, I. K. and Johnson, R. A. (2000). "A new family of
  power transformations to improve normality or symmetry,"
  Biometrika, 87, 954-959.
}

