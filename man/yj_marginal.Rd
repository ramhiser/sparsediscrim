\name{yj_marginal}
\alias{yj_marginal}
\title{Computes the marginal Yeo-Johnson (YJ) transformation to near-normality for
each class.}
\usage{
  yj_marginal(x, y)
}
\arguments{
  \item{x}{matrix of observations with observations on the
  rows and features on the columns.}

  \item{y}{vector of class labels for the observations
  (rows) in \code{x}.}
}
\value{
  named list, where each element is a list named by its
  class with the following elements: \itemize{ \item
  \code{lambda}: vector of fitted YJ transformation
  parameters for each feature vector.  \item \code{x}:
  matrix of the feature vectors transformed by the YJ
  transformation with the fitted parameter values. }
}
\description{
  For each class designated in \code{y}, we compute the
  marginal YJ transformation to near-normality. That is,
  within a given class, we estimate the YJ transformation
  parameter for each feature vector.
}
\details{
  The marginal YJ transformation to near-normality is a
  generalization of the marginal Box-Cox (BC)
  transformation that allows for observations to be
  unbounded on the real line, whereas the marginal BC
  transformation requires that observations be nonnegative.

  We use the \code{powerTransform} function from the
  \code{car} package to estimate the values for
  \code{lambda} and then use the \code{yjPower} function
  from the same package to perform the transformation.

  TODO: Provide the formulas for Yeo-Johnson transformation
  (given in paper) TODO: Cite Yeo-Johnson paper TODO:
  Reference the 'car' package to find the fitted 'lambda'
}
\examples{
yj_marginal(iris[, -5], iris$Species)
}

