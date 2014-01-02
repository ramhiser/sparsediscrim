\name{update_hdrda}
\alias{update_hdrda}
\title{Helper function to update tuning parameters for the HDRDA classifier}
\usage{
  update_hdrda(obj, lambda = 1, gamma = 0)
}
\arguments{
  \item{obj}{a \code{hdrda} object}

  \item{lambda}{a numeric value between 0 and 1,
  inclusively}

  \item{gamma}{a numeric value (nonnegative)}
}
\value{
  a \code{hdrda} object with updated estimates
}
\description{
  This function updates some of the quantities in the HDRDA
  classifier based on updated values of \code{lambda} and
  \code{gamma}. The update can greatly expedite
  cross-validation to examine a large grid of values for
  \code{lambda} and \code{gamma}.
}
\details{
  When \eqn{p < N}, the HDRDA covariance-matrix estimator
  is singular when \eqn{(lambda, gamma) = (0, 0)}. In this
  case, we set \eqn{gamma} to 0.01 to shrink the estimators
  slightly to ensure positive-definiteness.
}

