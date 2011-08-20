\name{predict_dlda}
\alias{predict_dlda}
\title{DLDA prediction of the class membership of a matrix of new observations.}
\usage{
  predict_dlda(obj, newdata)
}
\arguments{
  \item{obj}{trained DLDA object}

  \item{newdata}{matrix of observations to predict. Each
  row corresponds to a new observation.}
}
\value{
  list predicted class memberships of each row in newdata
}
\description{
  DLDA prediction of the class membership of a matrix of
  new observations.
}
\references{
  Dudoit, S., Fridlyand, J., & Speed, T. P. (2002).
  "Comparison of Discrimination Methods for the
  Classification of Tumors Using Gene Expression Data,"
  Journal of the American Statistical Association, 97, 457,
  77-87.
}

