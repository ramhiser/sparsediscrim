\name{tong_mean_shrinkage}
\alias{tong_mean_shrinkage}
\title{Tong et al. (2012)'s  Lindley-type Shrunken Mean Estimator}
\usage{
  tong_mean_shrinkage(x, r_opt = NULL)
}
\arguments{
  \item{x}{a matrix with \code{n} rows and \code{p}
  columns.}

  \item{r_opt}{the shrinkage coefficient. If \code{NULL}
  (default), we calculate the shrinkage coefficient with
  the formula given just above Equation 5 on page 533 and
  denoted by \eqn{\hat{r}_{opt}}. We allow the user to
  specify an alternative value to investigate better
  approximations.}
}
\value{
  vector of length \code{p} with the shrunken mean
  estimator
}
\description{
  TODO: Describe the shrunken mean TODO: Provide citation
  info.
}

