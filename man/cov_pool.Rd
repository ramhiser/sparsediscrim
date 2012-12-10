\name{cov_pool}
\alias{cov_pool}
\title{Computes the pooled maximum likelihood estimator for the common covariance
matrix from \eqn{K} multivariate normal populations.}
\usage{
  cov_pool(x, y, shrink = FALSE)
}
\arguments{
  \item{x}{data matrix with \code{n} observations and
  \code{p} feature vectors}

  \item{y}{class labels for observations (rows) in
  \code{x}}
}
\value{
  pooled sample covariance matrix of size \eqn{p \times p}
}
\description{
  TODO: Description
}
\examples{
cov_pool(iris[,-5], iris$Species)
}

