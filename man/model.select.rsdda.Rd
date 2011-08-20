\name{model.select.rsdda}
\alias{model.select.rsdda}
\title{RSDDA model selection.}
\usage{
  model.select.rsdda(object, grid_size = 5, k = 1)
}
\arguments{
  \item{object}{trained RSDDA object}

  \item{grid_size}{the grid size of lambda_grid considered
  in estimating the regularization parameter lambda}

  \item{k}{Leave-k-out crossvalidation is used to estimate
  the optimal value of lambda.}
}
\value{
  rsdda object with 'optimal' value of lambda
}
\description{
  Estimates the optimal value of lambda used in RSDDA by
  leave-k-out crossvalidation.
}
\details{
  of the pooled variances, we use a shrinkage method from
  Pang et al. (2009). This yields the SDQDA classifier.

  The RSDDA classifier is a regularized version of SDQDA,
  where the off-diagonal elements of the pooled sample
  covariance matrix are set to zero.
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

