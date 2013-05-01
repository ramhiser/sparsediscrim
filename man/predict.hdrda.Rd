\name{predict.hdrda}
\alias{predict.hdrda}
\title{Predicts the class membership of a matrix of unlabeled observations with a
trained HDRDA classifier.}
\usage{
  predict.hdrda(obj, newdata)
}
\arguments{
  \item{obj}{object of type \code{hdrda} that contains the
  trained HDRDA classifier}

  \item{newdata}{matrix containing the unlabeled
  observations to classify. Each row corresponds to a new
  observation.}
}
\value{
  list with predicted class and discriminant scores for
  each of the K classes
}
\description{
  For a given \code{hdrda} object, we predict the class of
  each observation (row) of the the matrix given in
  \code{newdata}.
}

