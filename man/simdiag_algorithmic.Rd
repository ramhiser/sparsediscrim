\name{simdiag_algorithmic}
\alias{simdiag_algorithmic}
\title{Simultaneously diagonalizes the sample covariance matrices for a data set
generated from \code{K} classes with one of several algorithmic approaches.}
\usage{
  simdiag_algorithmic(x, y,
    method = c("asfari", "jade", "jedi", "qdiag", "ffdiag", "jadiag", "uwedge"),
    tol = .Machine$double.eps, max_iter = 250)
}
\arguments{
  \item{x}{matrix containing the training data. The rows
  are the sample observations, and the columns are the
  features.}

  \item{y}{vector of class labels for each training
  observation}

  \item{method}{character. The SimDiag algorithm to employ.
  See Details for more information.}

  \item{tol}{a tolerance value below which covergence is
  reached for the specified algorithm.}
}
\value{
  a list containing: \itemize{ \item \code{Q}: simultaneous
  diagonalizing matrix of size \eqn{p \times p} \item
  \code{x}: the transformed data matrix, \code{x}.  \item
  \code{method}: string with the name of the employed
  SimDiag algorithm }
}
\description{
  This function determines the matrix, \code{Q}, that
  (nearly) simultaneously diagonalizes the sample
  covariance matrices for a data set generated from
  \code{K} classes.
}
\details{
  The returned simultaneous diagonalizing matrix, \code{Q},
  is constructed such that

  \deqn{Q \widehat{\\Sigma_k} Q' \approx D_k}

  for all \eqn{k = 1, \ldots, K}, where
  \eqn{\widehat{\\Sigma_k}} is the maximum likelihood
  estimator (under normality) for the \eqn{k}th class
  covariance matrix, \eqn{\\Sigma_k}, where \eqn{D_k} is a
  diagonal matrix.

  The returned simultaneously-diagonalizing matrix,
  \code{Q}, is of size \eqn{p \times p}.

  TODO: Briefly discuss each algorithmic approach. TODO:
  Provide reference to implementations in the 'jointDiag'
  package.
}
\examples{
# We use the well-known Iris data set to demonstrate four algorithmic
# approaches to SimDiag.
x <- iris[, -5]
y <- iris[, 5]

# Each of the methods return a list. We are interested in the SimDiag'd
# matrix, 'x'.
asfari_x <- simdiag_algorithmic(x = x, y = y, method = "asfari")$x
jedi_x <- simdiag_algorithmic(x = x, y = y, method = "jedi")$x
qdiag_x <- simdiag_algorithmic(x = x, y = y, method = "qdiag")$x
ffdiag_x <- simdiag_algorithmic(x = x, y = y, method = "ffdiag")$x

# Now, we should have that the covariance matrices for each class are nearly
# simultaneously diagonalized. Our goal is to verify that. To do so, we start
# off by computing each class's sample covariance matrix MLE for the
# SimDiag'd data. We also do the same for the original data.
orig_cov <- cov_list(x, y)
asfari_cov <- cov_list(asfari_x, y)
jedi_cov <- cov_list(jedi_x, y)
qdiag_cov <- cov_list(qdiag_x, y)
ffdiag_cov <- cov_list(ffdiag_x, y)

# To determine how well each SimDiag algorithm performed, we compute the
# Frobenius norm between each SimDiag'd covariance matrix and the
# corresponding matrix with only its diagonal entries. If the matrix were
# fully SimDiag'd, then the compute Frobenius norm should be equal to 0.
frob_diag <- function(x) {
  norm(x - diag(diag(x)), type = "F")
}

# We determine that the method that performed best is the one that has the
# minimum mean of the Frobenius norms. Note that this is a typical objective
# function for SimDiag methods. Note that the \\code{qdiag} method achieves
# the minimum value.
mean(sapply(orig_cov, frob_diag))
mean(sapply(asfari_cov, frob_diag))
mean(sapply(jedi_cov, frob_diag))
mean(sapply(qdiag_cov, frob_diag))
mean(sapply(ffdiag_cov, frob_diag))
}

