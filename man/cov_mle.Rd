\name{cov_mle}
\alias{cov_mle}
\title{Computes the maximum likelihood estimator for the sample covariance matrix
under the assumption of multivariate normality.}
\usage{
  cov_mle(x)
}
\arguments{
  \item{x}{data matrix with \code{n} observations and
  \code{p} feature vectors}
}
\value{
  sample covariance matrix of size \eqn{p \times p}
}
\description{
  For a sample matrix, \code{x}, we compute the sample
  covariance matrix of the data as the maximum likelihood
  estimator (MLE) of the population covariance matrix.
}

