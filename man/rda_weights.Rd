\name{rda_weights}
\alias{rda_weights}
\title{Computes the observation weights for each class for the RDA classifier}
\usage{
  rda_weights(x, y, lambda = 1)
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
  list containing the observations for each class given in
  \code{y}
}
\description{
  TODO: Add description
}

