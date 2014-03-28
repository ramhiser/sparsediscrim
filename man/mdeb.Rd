\name{mdeb}
\alias{mdeb}
\alias{mdeb.default}
\alias{mdeb.formula}
\alias{predict.mdeb}
\alias{print.mdeb}
\title{The Minimum Distance Empirical Bayesian Estimator (MDEB) classifier}
\usage{
  mdeb(x, ...)

  \method{mdeb}{default} (x, y, prior = NULL)

  \method{mdeb}{formula} (formula, data, prior = NULL)

  \method{print}{mdeb} (x, ...)

  \method{predict}{mdeb} (object, newdata)
}
\arguments{
  \item{x}{matrix containing the training data. The rows
  are the sample observations, and the columns are the
  features.}

  \item{y}{vector of class labels for each training
  observation}

  \item{prior}{vector with prior probabilities for each
  class. If NULL (default), then equal probabilities are
  used. See details.}

  \item{formula}{A formula of the form \code{groups ~ x1 +
  x2 + ...} That is, the response is the grouping factor
  and the right hand side specifies the (non-factor)
  discriminators.}

  \item{data}{data frame from which variables specified in
  \code{formula} are preferentially to be taken.}

  \item{x}{object to print}

  \item{...}{unused}

  \item{object}{trained mdeb object}

  \item{newdata}{matrix of observations to predict. Each
  row corresponds to a new observation.}

  \item{...}{additional arguments}
}
\value{
  \code{mdeb} object that contains the trained MDEB
  classifier

  list predicted class memberships of each row in newdata
}
\description{
  Given a set of training data, this function builds the
  MDEB classifier from Srivistava and Kubokawa (2007). The
  MDEB classifier is an adaptation of the linear
  discriminant analysis (LDA) classifier that is designed
  for small-sample, high-dimensional data. Rather than
  using the standard maximum likelihood estimator of the
  pooled covariance matrix, Srivastava and Kubokawa (2007)
  have proposed an Empirical Bayes estimator where the
  eigenvalues of the pooled sample covariance matrix are
  shrunken towards the identity matrix: the shrinkage
  constant has a closed form and is quick to calculate.

  Summarizes the trained mdeb classifier in a nice manner.

  The MDEB classifier from Srivistava and Kubokawa (2007)
  is an adaptation of the linear discriminant analysis
  (LDA) classifier that is designed for small-sample,
  high-dimensional data. Rather than using the standard
  maximum likelihood estimator of the pooled covariance
  matrix, Srivastava and Kubokawa (2007) have proposed an
  Empirical Bayes estimator where the eigenvalues of the
  pooled sample covariance matrix are shrunken towards the
  identity matrix: the shrinkage constant has a closed form
  and is quick to calculate
}
\details{
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
\examples{
n <- nrow(iris)
train <- sample(seq_len(n), n / 2)
mdeb_out <- mdeb(Species ~ ., data = iris[train, ])
predicted <- predict(mdeb_out, iris[-train, -5])$class

mdeb_out2 <- mdeb(x = iris[train, -5], y = iris[train, 5])
predicted2 <- predict(mdeb_out2, iris[-train, -5])$class
all.equal(predicted, predicted2)
}
\keyword{internal}

