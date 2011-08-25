\name{simdiag}
\alias{simdiag}
\title{Simultaneously diagonalizes two real symmetric matrices.}
\usage{
  simdiag(A, B, dim_reduce = F, k = NULL, eigen_pct = NULL,
  tol = sqrt(.Machine$double.eps))
}
\arguments{
  \item{A}{p x p real symmetric matrix}

  \item{B}{p x p real symmetric matrix}

  \item{dim_reduce}{logical flag: Should we reduce the
  dimension of the simultaneous diagonalization?}

  \item{k}{the reduced dimension of the simultaneous
  diagonalizer.}

  \item{eigen_pct}{The percentage of eigenvalues to
  preserve.}

  \item{tol}{a value indicating the magnitude below which
  eigenvalues are considered 0.}
}
\value{
  a list containing: the matrix Q that is the p x k
  simultaneous diagonalizing matrix the selected value of k
}
\description{
  This function simultaneously diagaonlizes two real
  symmetric matrices. If the matrices are singular, we
  employ a dimension reduction method that simultaneously
  diagonalizes them in a lower dimension.
}
\details{
  If the user selects the dimension reduction option, the
  user may specify a value of k. This is potentially useful
  if a scree plot is utilized or a manual threshold is
  desired. Otherwise, we select k automatically to be the
  number of non-zero eigenvalues. To determine what
  constitutes 'non-zero', we select all of the eigenvalues
  larger than the specified tolerance, 'tol'.

  The eigen_pct is the percentage of eigenvalues to use, by
  which we mean the cumulative sum of the eigenvalues
  (scaled by the sum of the eigenvalues). This is similar
  to the typical automatic approach to Principal Components
  Analysis.

  If the user does not wish to use the dimension reduction
  method, then k is automatically computed as the number of
  eigenvalues of B. Hence, in this case, k = p.

  The returned simultaneous diagonalizing matrix, Q, is
  constructed such that

  Q' A Q = D, and Q' B Q = I_k,

  where I_k denotes the k-dimensional identity matrix and
  Lambda is the diagonal matrix with diagonal entries as
  the characteristic roots (generalized eigenvalues) of
  det(B - lambda A) = 0.
}

