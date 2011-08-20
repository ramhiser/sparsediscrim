\name{predict_sdlda}
\alias{predict_sdlda}
\title{SDLDA prediction of the class membership of a matrix of new observations.}
\usage{
  predict_sdlda(obj, newdata)
}
\arguments{
  \item{obj}{trained SDLDA object}

  \item{newdata}{matrix of observations to predict. Each
  row corresponds to a new observation.}
}
\value{
  list predicted class memberships of each row in newdata
}
\description{
  of the pooled variances, we use a shrinkage method from
  Pang et al. (2009). This yields the SDLDA classifier.
}
\references{
  Dudoit, S., Fridlyand, J., & Speed, T. P. (2002).
  "Comparison of Discrimination Methods for the
  Classification of Tumors Using Gene Expression Data,"
  Journal of the American Statistical Association, 97, 457,
  77-87.

  Pang, H., Tong, T., & Zhao, H. (2009). "Shrinkage-based
  Diagonal Discriminant Analysis and Its Applications in
  High-Dimensional Data," Biometrics, 65, 4, 1021-1029.
}

