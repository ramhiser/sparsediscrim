\name{diag_estimates}
\alias{diag_estimates}
\title{Computes estimates and ancillary information for diagonal classifiers}
\usage{
  diag_estimates(x, y, prior = NULL, pool = FALSE)
}
\arguments{
  \item{matrix}{containing the training data. The rows are
  the sample observations, and the columns are the
  features.}

  \item{y}{vector of class labels for each training
  observation}

  \item{prior}{vector with prior probabilities for each
  class. If NULL (default), then equal probabilities are
  used. See details.}

  \item{pool}{logical value. If TRUE, calculates the pooled
  sample variances for each class.}
}
\value{
  list of MLEs and necessary ancillary information
}
\description{
  Computes the maximum likelihood estimators (MLEs) for
  each class under the assumption of multivariate normality
  for each class. Also, computes ancillary information
  necessary for classifier summary, such as sample size,
  the number of features, etc.
}
\details{
  This function computes the common estimates and ancillary
  information used in all of the diagonal classifiers in
  the \code{diagdiscrim} package.

  The matrix of training observations are given in
  \code{x}. The rows of \code{x} contain the sample
  observations, and the columns contain the features for
  each training observation.

  The vector of class labels given in \code{y} are coerced
  to a \code{factor}. The length of \code{y} should match
  the number of rows in \code{x}.

  An error is thrown if a given class has less than 2
  observations because the variance for each feature within
  a class cannot be estimated with less than 2
  observations.

  The vector, \code{prior}, contains the \emph{a priori}
  class membership for each class. If \code{prior} is NULL
  (default), the class membership probabilities are
  estimated as the sample proportion of observations
  belonging to each class. Otherwise, \code{prior} should
  be a vector with the same length as the number of classes
  in \code{y}. The \code{prior} probabilties should be
  nonnegative and sum to one.
}

