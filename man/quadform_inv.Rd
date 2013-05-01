\name{quadform_inv}
\alias{quadform_inv}
\title{Quadratic Form of the inverse of a matrix and a vector}
\usage{
  quadform_inv(A, x)
}
\arguments{
  \item{A}{matrix that is p x p and nonsingular}

  \item{x}{vector of length p}
}
\value{
  quadratic form x' A^{-1} x
}
\description{
  We compute the quadratic form of a vector and the inverse
  of a matrix in an efficient manner. Let x be a real
  vector of length p, and let A be a p x p nonsingular
  matrix. Then, we compute the quadratic form q by
}
\details{
  q = x' A^{-1} x.

  A naive way to compute the quadratic form is to
  explicitly write

  t(x) %*% solve(A) %*% x,

  but for large p, this operation is inefficient. We
  provide a more efficient method below.

  Note that we have adapted the code from:
  \url{http://tolstoy.newcastle.edu.au/R/help/05/11/14989.html}
}

