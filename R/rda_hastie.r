#' Regularized Discriminant Analysis (RDA) as defined by Hastie, Tibshirani, and
#' Friedman (2008).
#'
#' Given a set of training data, this function builds the RDA classifier from
#' Hastie et al. (2008) in the Elements of Statistical Learning textbook
#'
#' The matrix of training observations are given in \code{x}. The rows of \code{x}
#' contain the sample observations, and the columns contain the features for each
#' training observation.
#'
#' The vector of class labels given in \code{y} are coerced to a \code{factor}.
#' The length of \code{y} should match the number of rows in \code{x}.
#'
#' The vector, \code{prior}, contains the \emph{a priori} class membership for
#' each class. If \code{prior} is NULL (default), the class membership
#' probabilities are estimated as the sample proportion of observations belonging
#' to each class. Otherwise, \code{prior} should be a vector with the same length
#' as the number of classes in \code{y}. The \code{prior} probabilties should be
#' nonnegative and sum to one. The order of the prior probabilties is assumed to
#' match the levels of factor(y).
#'
#' TODO: Discuss the shrinkage
#' TODO: Add Fast SVD trick
#'
#' @export
#' @references Friedman, J. H. (1989), "Regularized Discriminant Analysis,"
#' Journal of American Statistical Association, 84, 405, 165-175.
#' \url{http://www.jstor.org/pss/2289860} (Requires full-text access).
#' @param x matrix containing the training data. The rows are the sample
#' observations, and the columns are the features.
#' @param y vector of class labels for each training observation
#' @param gamma the RDA shrinkage parameter. Must be between 0 and 1,
#' inclusively.
#' @param prior vector with prior probabilities for each class. If \code{NULL}
#' (default), then the sample proportion of observations belonging to each class
#' equal probabilities are used. See details.
#' @param tol a threshold for determining 'near-zero' eigenvalues.
#' @return \code{rda_hastie} object that contains the trained RDA classifier
rda_hastie <- function(x, ...) {
  UseMethod("rda_hastie")
}

#' @rdname rda_hastie
#' @method rda_hastie default
#' @S3method rda_hastie default
rda_hastie.default <- function(x, y, gamma = 0, prior = NULL, tol = 1e-6) {
  x <- as.matrix(x)
  y <- as.factor(y)
  gamma <- as.numeric(gamma)
  shrinkage <- match.arg(shrinkage)

  if (gamma < 0 || gamma > 1) {
    stop("The value for 'gamma' must be between 0 and 1, inclusively.")
  }

  obj <- regdiscrim_estimates(x = x, y = y, prior = prior)

  rda_cov <- gamma * cov_pool
  diag(rda_cov) <- diag(rda_cov) + (1 - gamma) * diag(cov_pool)
  rda_cov_eigen <- eigen(rda_cov, symmetric = TRUE)

  # Replaces any eigenvalues less than 'tol' with 'tol'. This ensures that the
  # covariance-matrix estimator is nonsingular.
  rda_cov_eigen$values <- with(rda_cov_eigen, replace(values, values < tol, tol))

  obj$est <- lapply(obj$est, function(class_est) {
    class_est$cov <- rda_cov
    class_est$cov_eigen <- rda_cov_eigen

    class_est
  })

  # Creates an object of type 'rda_hastie' and adds the 'match.call' to the object
  obj$call <- match.call()
  class(obj) <- "rda_hastie"
	
  obj
}

#' @param formula A formula of the form \code{groups ~ x1 + x2 + ...} That is,
#' the response is the grouping factor and the right hand side specifies the
#' (non-factor) discriminators.
#' @param data data frame from which variables specified in \code{formula} are
#' preferentially to be taken.
#' @rdname rda_hastie
#' @method rda_hastie formula
#' @S3method rda_hastie formula
rda_hastie.formula <- function(formula, data, gamma = 0, prior = NULL,
                               tol = 1e-6) {
  # The formula interface includes an intercept. If the user includes the
  # intercept in the model, it should be removed. Otherwise, errors and doom
  # happen.
  # To remove the intercept, we update the formula, like so:
  # (NOTE: The terms must be collected in case the dot (.) notation is used)
  formula <- no_intercept(formula, data)
  
  mf <- model.frame(formula = formula, data = data)
  x <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)

  est <- rda_hastie.default(x = x, y = y, gamma = gamma, prior = prior,
                            tol = tol)
  est$call <- match.call()
  est$formula <- formula
  est
}

#' Predicts the class membership of a matrix of unlabeled observations with a
#' trained RDA classifier.
#'
#' For a given \code{rda_hastie} object, we predict the class of each observation (row)
#' of the the matrix given in \code{newdata}.
#'
#' @export
#' @param obj object of type \code{rda_hastie} that contains the trained RDA classifier
#' @param newdata matrix containing the unlabeled observations to classify. Each
#' row corresponds to a new observation.
#' @return list with predicted class and discriminant scores for each of the K
#' classes
predict.rda_hastie <- function(obj, newdata) {
  if(is.vector(newdata)) {
    newdata <- matrix(newdata, nrow = 1)
  }
  newdata <- as.matrix(newdata)
  
  scores <- sapply(obj$est, function(class_est) {
    cov_inv <- with(class_est$cov_eigen,
                    vectors %*% tcrossprod(diag(1 / values), vectors))
    apply(newdata, 1, function(obs) {
      x_centered <- as.vector(obs - class_est$xbar)
      obs_quadform <- quadform(A = cov_inv, x = x_centered)
      score <- obs_quadform - 2 * log(class_est$prior)
      score
    })
  })

  if(is.vector(scores)) {
    min_scores <- which.min(scores)
    posterior <- exp(scores) / (1 + sum(exp(scores)))
  } else {
    min_scores <- apply(scores, 1, which.min)
    posterior <- t(apply(scores, 1, function(score_i) {
      exp(score_i) / (1 + sum(exp(score_i)))
    }))
  }

  class <- with(obj, factor(groups[min_scores], levels = groups))

  list(class = class, scores = scores, posterior = posterior)
}
