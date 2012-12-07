\name{is_diagonal}
\alias{is_diagonal}
\title{Tests if a matrix is diagonal.}
\usage{
  is_diagonal(x, tol = 1e-08)
}
\arguments{
  \item{x}{matrix}

  \item{tol}{tolerance value}
}
\value{
  logical value indicating whether the matrix \code{x} is
  diagonal
}
\description{
  This function tests if the matrix \code{x} is diagonal by
  computing the Frobenius norm of the difference between
  \code{x} and the matrix consisting of the diagonal
  elements of \code{x} is less than some tolerance.
}

