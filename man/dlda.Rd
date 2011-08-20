\name{dlda}
\alias{dlda}
\title{Diagonalized Linear Discriminant Analysis (DLDA)}
\usage{
  dlda(x, y)
}
\arguments{
  \item{x}{training data in matrix form.}

  \item{y}{labels of the training data.}
}
\value{
  dlda obj
}
\description{
  Given a set of training data, this function builds the
  DLDA classifier, which is often attributed to Dudoit et
  al. (2002).
}
\references{
  Dudoit, S., Fridlyand, J., & Speed, T. P. (2002).
  "Comparison of Discrimination Methods for the
  Classification of Tumors Using Gene Expression Data,"
  Journal of the American Statistical Association, 97, 457,
  77-87.
}

