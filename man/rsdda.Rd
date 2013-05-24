\name{rsdda}
\alias{rsdda}
\title{Regularized Shrinkage-based Diagonal Discriminant Analysis (RSDDA)}
\usage{
  rsdda(x, y, num_alphas = 101)
}
\arguments{
  \item{x}{training data in matrix form.}

  \item{y}{labels of the training data.}

  \item{num_alphas}{The number of values used to find the
  optimal amount of shrinkage.}
}
\value{
  rsdda obj
}
\description{
  Given a set of training data, this function builds the
  RSDDA classifier. To improve the estimation of the pooled
  variances, we use a shrinkage method from Pang et al.
  (2009).  This yields the SDQDA classifier.
}
\details{
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

