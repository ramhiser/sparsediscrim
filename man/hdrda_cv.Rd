\name{hdrda_cv}
\alias{hdrda_cv}
\title{Helper function to optimize the HDRDA classifier via cross-validation}
\usage{
  hdrda_cv(x, y, num_folds = 10, num_lambda = 21,
    num_gamma = 7, shrinkage_type = c("ridge", "convex"),
    ...)
}
\arguments{
  \item{x}{matrix containing the training data. The rows
  are the sample observations, and the columns are the
  features.}

  \item{y}{vector of class labels for each training
  observation}

  \item{num_folds}{the number of cross-validation folds.}

  \item{num_lambda}{The number of values of \code{lambda}
  to consider}

  \item{num_gamma}{The number of values of \code{gamma} to
  consider}

  \item{shrinkage_type}{the type of covariance-matrix
  shrinkage to apply. By default, a ridge-like shrinkage is
  applied. If \code{convex} is given, then shrinkage
  similar to Friedman (1989) is applied. See Ramey et al.
  (2013) for details.}

  \item{...}{Additional arguments passed to
  \code{\link{hdrda}}.}
}
\value{
  list containing the HDRDA model that minimizes
  cross-validation as well as a \code{data.frame} that
  summarizes the cross-validation results.
}
\description{
  For a given data set, we apply cross-validation (cv) to
  select the optimal HDRDA tuning parameters.
}
\details{
  The number of cross-validation folds is given in
  \code{num_folds}.
}

