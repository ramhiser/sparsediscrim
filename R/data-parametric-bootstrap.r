#' Generates data from \code{K} multivariate data populations via a parametric
#' bootsrapping technique with covariance matrix shrinkage.
#'
#' This function generates \code{K} multivariate data sets, where each data set
#' is generated with a mean vector and a covariance matrix calculated from the
#' data provided in \code{x}. The data are returned as a list containing a
#' single matrix \code{x} along with a vector of class labels \code{y} that
#' indicates class membership.
#' 
#' The covariance matrices for each class are constructed using the maximum
#' likelihood estimators under multivariate normality. For high-dimensional
#' data, these covariance matrices are likely singular, in which case we appy
#' shrinkage to the covariance matrix estimators.
#'
#' We use the \code{\link{car}} package's implementation of the Box-Cox and
#' Yeo-Johnson transformation methods. To compute the pseudo-likelihood
#' estimators, we use the \code{\link[car]{powerTransform}} function. In this
#' function, the author uses the \code{\link{optim}} function to numerically
#' optimize the pseudo-likelihood functions for the given data. By default, the
#' lower and upper bounds are \code{-Inf} and \code{Inf}, respectively, which
#' can yield numerically unstable estimates in the optimization
#' search. Practically, we wish to consider only values between -4 and 4, but we
#' allow the user to alter these values by way of the \code{optim_lower} and
#' \code{optim_upper} arguments, respectively. See
#' \code{\link[car]{powerTransform}} for more details.
#'
#' @export
#' @param n vector of the sample sizes of each class to generate. The length of
#' \code{n} should match the number of classes given in \code{y}.
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
#' @param optim_method the specified numerical method to use for the
#' transformation parameter estimation. By default, we use the \code{Nelder-Mead}
#' option; for other values, see the \code{method} argument for the \code{optim}
#' function.
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
#' TODO: Add examples
#' TODO: Import functions from 'mvtnorm' and 'car' packages
boot_parametric <- function(n, x, y, gamma = 1,
                            transformation = c("none", "Box-Cox", "Yeo-Johnson")) {
  # TODO: Update 'require' statements to @import
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
    stop("The Box-Cox transformation has not yet been implemented. Use
          Yeo-Johnson instead.")
  }

  # If the Yeo-Johnson transformation is selected, we marginally estimate the
  # power parameters 'lambda' for each class.
  if (transformation == "Yeo-Johnson") {
    yj_out <- yj_marginal(x = x, y = y)
    
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
  
  # If the Yeo-Johnson transformation is selected, we now apply the inverse of
  # the Yeo-Johnson transformation based on the estimated values of 'lambda' for
  # each class. These values are stored in the list 'yj_lambda'.
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

#' Computes the marginal Yeo-Johnson transformation to near-normality for each
#' class.
#'
#' For each class designated in \code{y}, we compute the marginal Yeo-Johnson
#' (YJ) transformation to near-normality. That is, within a given class, we
#' estimate the YJ transformation parameter for each feature vector.
#'
#' The marginal YJ transformation to near-normality is a generalization of the
#' marginal Box-Cox (BC) transformation that allows for observations to be
#' unbounded on the real line, whereas the marginal BC transformation requires
#' that observations be nonnegative.
#'
#' We use the \code{\link[car]{powerTransform}} function to estimate the values
#' for \code{lambda} and then use the \code{\link[car]{yjPower}} function to
#' perform the transformation.
#'
#' @export
#' @param x matrix of observations with observations on the rows and features on
#' the columns.
#' @param y vector of class labels for the observations (rows) in \code{x}.
#' @return named list, where each element is a list named by its class with the
#' following elements:
#' \itemize{
#'   \item \code{lambda}: vector of fitted YJ transformation parameters for each
#' feature vector.
#'   \item \code{x}: matrix of the feature vectors transformed by the YJ
#' transformation with the fitted parameter values.
#' }
#'
#' @references
#' Yeo, I. K. and Johnson, R. A. (2000). "A new family of power transformations
#' to improve normality or symmetry," Biometrika, 87, 954-959.
#' 
#' @examples
#' yj_marginal(iris[, -5], iris$Species)
yj_marginal <- function(x, y) {
  # TODO: Update 'require' statements to @import
  require('car')
  tapply(seq_along(y), y, function(i) {
    x_i <- x[i, ]
    n_i <- nrow(x_i)
    lambda <- apply(x_i, 2, function(xi_col) {
      estimateTransform(X = rep(1, n_i), Y = xi_col, family = "yjPower")$lambda
    })
    list(
         lambda = lambda,
         x = unname(as.matrix(yjPower(x_i, lambda = lambda)))
        )
  })
}


#' Computes the inverse of Yeo-Johnson transformed random variate(s)
#'
#' For a random variate that has been transformed with the Yeo-Johnson
#' transformation to near-normality, we invert the random variate back to its
#' original scale in this function. The inversion requires the transformation
#' value \code{lambda} applied initially.
#'
#' The YJ transformation to near-normality is a generalization of the Box-Cox
#' (BC) transformation that allows for observations to be unbounded on the real
#' line, whereas the BC transformation requires that observations be
#' nonnegative.
#'
#' The value for \code{lambda} is typically obtained using the
#' \code{\link[car]{powerTransform}} function. The transformed random variate is
#' usually obtained from the \code{\link[car]{yjPower}} function.
#'
#' @export
#' @param y vector of random variate(s)
#' @param lambda numeric Yeo-Johnson power transformation parameter value
#' @param tol numeric tolerance value that determines numerical equality
#' @return the inverse of the Yeo-Johnson random variate, \code{y}
#'
#' @references
#' Yeo, I. K. and Johnson, R. A. (2000). "A new family of power transformations
#' to improve normality or symmetry," Biometrika, 87, 954-959.
yj_inverse <- Vectorize(function(y, lambda, tol = sqrt(.Machine$double.eps)) {
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
}, vectorize.args = c("y", "lambda"))
