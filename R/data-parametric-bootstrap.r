#' Generates data from \code{K} multivariate data populations via a parametric
#' bootsrapping technique with covariance matrix shrinkage.
#'
#' This function generates \code{K} multivariate data sets, where each
#' data set is generated with a mean vector and a covariance matrix calculated
#' from the data provided in \code{x}. The data are returned as a single matrix,
#' \code{x}, along with a vector of class labels, \code{y}, that indicates class
#' membership.
#' 
#' TODO: Define the covariance matrix, \eqn{\Sigma_k}.
#' TODO: Define the shrinkage covariance matrix estimator.
#'
#' We use the \code{car} package's implementation of the Box-Cox and Yeo-Johnson
#' transformation methods. To compute the pseudo-likelihood estimators, we use
#' the \code{powerTransform} function in the \code{car} package. In this
#' function, the author uses the \code{optim} function to numerically optimize
#' the pseudo-likelihood functions for the given data. By default, the lower
#' and upper bounds are \code{-Inf} and \code{Inf}, respectively, which can
#' yield numerically unstable estimates in the optimization search. Practically,
#' we wish to only consider values between -5 and 5, but we allow the user to
#' alter these values by way of the \code{optim_lower} and \code{optim_upper}
#' arguments, respectively. See the \code{powerTransform} function in the
#' \code{car} package for more details.
#'
#' @export
#' @param n vector of the sample sizes of each class to generate. The
#' \code{length} of \code{n} should match the number of classes given in
#' \code{y}.
#' @param x matrix of observations with observations on the rows and features on
#' the columns
#' @param y vector of class labels for the observations (rows) in \code{x}.
#' @param gamma numeric value between 0 and 1, inclusively. The value shrinks
#' the sample covariance matrix for each class towards its diagonal with a weight
#' given by \code{gamma}.
#' @param transformation character. If \code{none} (default), no transformation is
#' applied to data in \code{x} before sampling from the Multivariate Normal (MVN)
#' distribution. If \code{Yeo-Johnson}, we first apply the Yeo-Johnson
#' near-normality transformation to each column (that is, marginally) before
#' generating data from the MVN distribution. After the data has been generated,
#' we apply the inverse transformation to the generated data to attempt to
#' restore the shape and scale of the underlying data.
#' @param optim_lower the lower bound for the values considered in the numerical
#' optimization function, \code{optim}, that is used to determine the
#' pseudo-likelihood transformation estimators. Ignored if \code{transformation}
#' is \code{none}.
#' @param optim_upper the lower bound for the values considered in the numerical
#' optimization function, \code{optim}, that is used to determine the
#' pseudo-likelihood transformation estimators. Ignored if \code{transformation}
#' is \code{none}.
#' @return named list with elements:
#' \itemize{
#'   \item \code{x}: matrix of observations with \code{sum(n)} rows and \code{p}
#' columns
#'   \item \code{y}: vector of class labels that indicates class membership for
#' each observation (row) in \code{x}.
#' }
#' @examples
#' TODO
boot_parametric <- function(n, x, y, gamma = 1,
                            transformation = c("none", "Box-Cox", "Yeo-Johnson"),
                            optim_lower = -5, optim_upper = 5) {
  require('mvtnorm')
  n <- as.integer(n)
  x <- as.matrix(x)
  y <- as.factor(y)
  transformation <- match.arg(transformation)

  if (gamma < 0 || gamma > 1) {
    stop("The value of 'gamma' must be between 0 and 1, inclusively.")
  }
  if (length(y) != nrow(x)) {
    stop("The length of 'y' must equal the number of rows of 'x'.")
  }
  if (length(n) != nlevels(y)) {
    stop("The length of 'n' must equal the number of classes in 'y'.")
  }

  if (transformation == "Box-Cox") {
    stop("The Box-Cox transformation has not yet been implemented. Use Yeo-Johnson instead.")
  }

  # If the Yeo-Johnson transformation is selected, we marginally estimate
  # the power parameters 'lambda' for each class. We use the 'powerTransform'
  # function from the 'car' package to estimate the values for 'lambda' and
  # then use the 'yjPower' function from the same package to perform the
  # transformation.
  if (transformation == "Yeo-Johnson") {
    yj_out <- tapply(seq_along(y), y, function(i) {
      lambda <- powerTransform(x[i, ], family = "yjPower", lower = optim_lower,
                               upper = optim_upper)$lambda
      list(
           lambda = lambda,
           x = as.matrix(yjPower(x[i, ], lambda = lambda))
      )
    })
    # The list 'yj_lambda' stores the transformation parameters for each class.
    yj_lambda <- lapply(yj_out, function(z) z$lambda)

    # The matrix 'x' contains the transformed data.
    x <- do.call(rbind, lapply(yj_out, function(z) z$x))
  }

  xbars <- tapply(seq_along(y), y, function(i) {
    colMeans(x[i, ])
  })
  covs_shrunken <- tapply(seq_along(y), y, function(i) {
    cov_shrink_diag(x[i, ], gamma = gamma)
  })
  sample_sizes <- as.list(n)
  
  x <- mapply(function(n, mu, sigma) {
      rmvnorm(n = n, mean = mu, sigma = sigma)
    }, sample_sizes, xbars, covs_shrunken, SIMPLIFY = FALSE)
  x <- unname(do.call(rbind, x))
      
  y <- mapply(function(n, label) {
    rep.int(label, n)
  }, sample_sizes, levels(y), SIMPLIFY = FALSE, USE.NAMES = FALSE)
  y <- as.factor(do.call(c, y))
  
  # If the Yeo-Johnson transformation is selected, we now apply the inverse
  # of the Yeo-Johnson transformation based on the estimated values of 'lambda'
  # for each class. These values are stored in the list 'yj_lambda'.
  if (transformation == "Yeo-Johnson") {
    # The row indices for each class.
    class_indices <- tapply(seq_along(y), y, identity)

    # For each class, we apply the inverse of the Yeo-Johnson transformation to
    # the generated data.
    x <- mapply(function(i, lambda) {
      apply(x[i, ], 1, yj_inverse, lambda = lambda)
    }, class_indices, yj_lambda, SIMPLIFY = FALSE)
    x <- t(do.call(cbind, x))
  }

  list(x = x, y = y)
}

#' Computes the inverse of Yeo-Johnson transformed random variate(s)
#'
#' TODO: Briefly describe Yeo-Johnson is generalization of Box-Cox
#' TODO: Provide the formulas for Yeo-Johnson transformation (given in paper)
#' TODO: Provide the inverse formulas for Yeo-Johnson (on whiteboard and iPhone)
#' TODO: Cite Yeo-Johnson paper
#' TODO: Reference the 'car' package to find the fitted 'lambda'
#'
#' @export
#' @param y vector of random variate(s)
#' @param lambda numeric Yeo-Johnson power transformation parameter value
#' @param tol numeric tolerance value that determines numerical equality
#' @return the inverse of the Yeo-Johnson random variate, \code{y}
yj_inverse <- function(y, lambda, tol = sqrt(.Machine$double.eps)) {
  if (y >= 0) {
    # Case: lambda = 0 (numerically)
    if (abs(lambda) < tol) {
      x <- exp(y) - 1
    } else { # Case: Lambda != 0
      x <- (lambda * y + 1)^(1 / lambda) - 1
    }
  } else { # y is negative
    # Case: lambda = 2 (numerically)
    if (abs(lambda - 2) < tol) {
      x <- 1 - exp(-y)
    } else { # Case: Lambda != 2
      x <- -(-(2 - lambda) * y + 1)^(1 / (2 - lambda)) + 1
    }
  }
  x
}
yj_inverse <- Vectorize(yj_inverse, vectorize.args = c("y", "lambda"))
