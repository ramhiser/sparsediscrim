\name{predict.dlda}
\alias{predict.dlda}
\title{DLDA prediction of the class membership of a matrix of new observations.}
\usage{
  predict.dlda(object, newdata)
}
\arguments{
  \item{object}{trained DLDA object}

  \item{newdata}{matrix of observations to predict. Each
  row corresponds to a new observation.}
}
\value{
  list predicted class memberships of each row in newdata
}
\description{
  The DLDA classifier is a modification to LDA, where the
  off-diagonal elements of the pooled sample covariance
  matrix are set to zero.
}
\references{
  Dudoit, S., Fridlyand, J., & Speed, T. P. (2002).
  "Comparison of Discrimination Methods for the
  Classification of Tumors Using Gene Expression Data,"
  Journal of the American Statistical Association, 97, 457,
  77-87.
}

