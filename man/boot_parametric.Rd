\name{boot_parametric}
\alias{boot_parametric}
\title{Generates data from \code{K} multivariate data populations via a parametric
bootsrapping technique with covariance matrix shrinkage.}
\usage{
  boot_parametric(n, x, y, gamma = 1,
    transformation = c("none", "Yeo-Johnson"))
}
\arguments{
  \item{n}{vector of the sample sizes of each class to
  generate. The \code{length} of \code{n} should match the
  number of classes given in \code{y}.}

  \item{x}{matrix of observations with observations on the
  rows and features on the columns}

  \item{y}{vector of class labels for the observations
  (rows) in \code{x}.}

  \item{gamma}{numeric value between 0 and 1, inclusively.
  The value shrinks the sample covariance matrix for each
  class towards its diagonal with a weight given by
  \code{gamma}.}

  \item{transformation}{character. If \code{none}
  (default), no transformation is applied to data in
  \code{x} before sampling from the Multivariate Normal
  (MVN) distribution. If \code{Yeo-Johnson}, we first apply
  the Yeo-Johnson near-normality transformation to each
  column (that is, marginally) before generating data from
  the MVN distribution. After the data has been generated,
  we apply the inverse transformation to the generated data
  to attempt to restore the shape and scale of the
  underlying data.}
}
\value{
  named list with elements: \itemize{ \item \code{x}:
  matrix of observations with \code{sum(n)} rows and
  \code{p} columns \item \code{y}: vector of class labels
  that indicates class membership for each observation
  (row) in \code{x}. }
}
\description{
  This function generates \code{K} multivariate data sets,
  where each data set is generated with a mean vector and a
  covariance matrix calculated from the data provided in
  \code{x}. The data are returned as a single matrix,
  \code{x}, along with a vector of class labels, \code{y},
  that indicates class membership.
}
\details{
  TODO: Define the covariance matrix, \eqn{\Sigma_k}. TODO:
  Define the shrinkage covariance matrix estimator.
}
\examples{
TODO
}

