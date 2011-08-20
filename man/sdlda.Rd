\name{sdlda}
\alias{sdlda}
\title{Shrinkage-based Diagonal Linear Discriminant Analysis (SDLDA)}
\usage{
  sdlda(x, y, num_alphas = 101)
}
\arguments{
  \item{x}{training data in matrix form.}

  \item{y}{labels of the training data.}

  \item{num_alphas}{The number of values used to find the
  optimal amount of shrinkage.}
}
\value{
  sdlda obj
}
\description{
  Given a set of training data, this function builds the
  DLDA classifier, which is often attributed to Dudoit et
  al. (2002). To improve the estimation of the pooled
  variances, we use a shrinkage method from Pang et al.
  (2009). This yields the SDLDA classifier.
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

