\name{cov_list}
\alias{cov_list}
\title{Computes the covariance matrix maximum likelihood estimators for each class
and returns a list.}
\usage{
  cov_list(x, y, shrink = FALSE)
}
\arguments{
  \item{x}{data matrix with \code{n} observations and
  \code{p} feature vectors}

  \item{y}{class labels for observations (rows) in
  \code{x}}

  \item{shrink}{logical. If \code{TRUE}, each covariance
  matrix is shrunken toward the a scaled identity matrix to
  obtain a ridge-like covariance matrix estimator. We use
  the MDEB estimator to determine the shrinkage term. By
  default, the argument is \code{FALSE}, and no shrinkage
  is applied.}
}
\value{
  list of the sample covariance matrices of size \eqn{p
  \times p} for each class given in \code{y}.
}
\description{
  For a sample matrix, \code{x}, we compute the sample
  covariance matrix for each class given in the vector,
  \code{y}.
}

