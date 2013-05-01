\name{rda_cov}
\alias{rda_cov}
\title{Calculates the RDA covariance-matrix estimators for each class}
\usage{
  rda_cov(x, y, lambda = 1)
}
\arguments{
  \item{x}{matrix containing the training data. The rows
  are the sample observations, and the columns are the
  features.}

  \item{y}{vector of class labels for each training
  observation}

  \item{lambda}{the RDA pooling parameter. Must be between
  0 and 1, inclusively.}
}
\value{
  list containing the RDA covariance-matrix estimators for
  each class given in \code{y}
}
\description{
  TODO: Add description
}

