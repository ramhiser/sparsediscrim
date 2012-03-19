\name{diag_shrinkage}
\alias{diag_shrinkage}
\title{Applies covariance matrix shrinkage to each class for the DLDA and DQDA
classifiers.}
\usage{
  diag_shrinkage(obj, shrinkage, pool = FALSE)
}
\arguments{
  \item{obj}{object of type \code{dlda} or \code{dqda}}

  \item{shrinkage}{list of length \eqn{K}. The \eqn{k}th
  list element is a diagonal matrix towards which the
  \eqn{k}th sample covariance matrix is shrunken.}

  \item{pool}{logical value. Add the shrinkage matrix to
  the diagonal pooled sample covariance matrix.}
}
\value{
  list of length \code{K} with MDEB shrinkage matrices
}
\description{
  For each class covariance matrix (which is assumed to be
  diagonal), we shrink the diagonal matrix towards a
  diagonal matrix, \eqn{D_k}. We assume that each sample
  covariance matrix, \eqn{\widehat{\Sigma}_k}, is diagonal
  for \eqn{k = 1, \ldots, K}, where \eqn{K} is the number
  of classes. We wish to shrink each covariance matrix
  towards a diagonal matrix, \eqn{D_k}. The matrix
  \eqn{D_k} is given as the \eqn{k}th list element in the
  list, given as an argument, \code{shrinkage}.
}

