#' Regularized Discriminant Analysis (RDA)
#'
#' Given a set of training data, this function builds the RDA classifier from
#' Friedman (1989). The RDA classifier generalized the Linear Discriminant
#' Analysis (LDA) and Quadratic Discriminant Analysis (QDA) classifiers by
#' introducing a biased covariance-matrix estimator with two tuning parameters
#' that control the amount of pooling and the amount of shrinkage.
#'
#' For large data sets, we recommend using the singular value decomposition (SVD)
#' "trick" to significantly improve the runtime reduce the dimension of the data
#' for estimation and model selection.
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
#' By default, we employ the original shrinkage factor from Friedman (1989) and
#' use the average of the eigenvalues. Alternatively, we allow a shrinkage factor
#' that is approximately equal to the nonzero eigenvalues from Srivastava and
#' Kubokawa (2007).
#' 
#' While training the RDA classifier, we compute the inverse of the RDA
#' covariance matrix estimator. For singular matrices, the inverse will result in
#' an error. In this case, we replace the 'near-zero' eigenvalues of the
#' covariance matrix with non-zero values. In particular, we replace the
#' eigenvalues below the tolerance level, \code{tol}, with the tolerance level.
#' For example, the default value of \code{tol} is 0.001. If we encounter an
#' eigenvalue that is near zero, we will replace it with 0.001 if the covariance
#' matrix estimator is singular.
#'
#' As discussed above, the covariance matrix estimator for each class can be
#' singular, which implies that its determinant is zero. In this case, we also
#' apply the aforementioned eigenvalue adjustment before computing the
#' determinant.
#'
#' @export
#' @references Friedman, J. H. (1989), "Regularized Discriminant Analysis,"
#' Journal of American Statistical Association, 84, 405, 165-175.
#' \url{http://www.jstor.org/pss/2289860} (Requires full-text access).
#' @param x matrix containing the training data. The rows are the sample
#' observations, and the columns are the features.
#' @param y vector of class labels for each training observation
#' @param lambda the RDA pooling parameter. Must be between 0 and 1, inclusively.
#' @param gamma the RDA shrinkage parameter. Must be between 0 and 1,
#' inclusively.
#' @param shrinkage What shrinkage factor should we employ? See details.
#' @param prior vector with prior probabilities for each class. If \code{NULL}
#' (default), then the sample proportion of observations belonging to each class
#' equal probabilities are used. See details.
#' @param tol a threshold for determining 'near-zero' eigenvalues.
#' @return \code{rda} object that contains the trained RDA classifier
rda <- function(x, ...) {
  UseMethod("rda")
}

#' @rdname rda
#' @method rda default
#' @S3method rda default
rda.default <- function(x, y, lambda = 1, gamma = 0,
                        shrinkage = c("Friedman", "Srivastava"), prior = NULL,
                        tol = 1e-6) {
  x <- as.matrix(x)
  y <- as.factor(y)
  lambda <- as.numeric(lambda)
  gamma <- as.numeric(gamma)
  shrinkage <- match.arg(shrinkage)

  if (lambda < 0 || lambda > 1) {
    stop("The value for 'lambda' must be between 0 and 1, inclusively.")
  }
  if (gamma < 0 || gamma > 1) {
    stop("The value for 'gamma' must be between 0 and 1, inclusively.")
  }

  obj <- regdiscrim_estimates(x = x, y = y, prior = prior)

  obj$est <- lapply(obj$est, function(class_est) {
    class_est$cov <- cov_rda(cov_k = class_est$cov, cov_pool = obj$cov_pool,
                             n_k = class_est$n, N = obj$N, lambda = lambda,
                             gamma = gamma, shrinkage = shrinkage)
    class_est$cov_eigen <- eigen(class_est$cov, symmetric = TRUE)
    # Replaces any eigenvalues less than 'tol' with 'tol'. This ensures that the
    # covariance-matrix estimator is nonsingular.
    class_est$cov_eigen$values <- with(class_est$cov_eigen,
                                       replace(values, values < tol, tol))

    class_est
  })

  # Creates an object of type 'rda' and adds the 'match.call' to the object
  obj$call <- match.call()
  class(obj) <- "rda"
	
  obj
}

#' @param formula A formula of the form \code{groups ~ x1 + x2 + ...} That is,
#' the response is the grouping factor and the right hand side specifies the
#' (non-factor) discriminators.
#' @param data data frame from which variables specified in \code{formula} are
#' preferentially to be taken.
#' @rdname rda
#' @method rda formula
#' @S3method rda formula
rda.formula <- function(formula, data, lambda = 1, gamma = 0,
                        shrinkage = c("Friedman", "Srivastava"), prior = NULL,
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

  est <- rda.default(x = x, y = y, lambda = lambda, gamma = gamma, prior = prior,
                     tol = tol)
  est$call <- match.call()
  est$formula <- formula
  est
}

#' Predicts the class membership of a matrix of unlabeled observations with a
#' trained RDA classifier.
#'
#' For a given \code{rda} object, we predict the class of each observation (row)
#' of the the matrix given in \code{newdata}.
#'
#' While training the RDA classifier, we compute the inverse of the RDA
#' covariance matrix estimator. For singular matrices, the inverse will result in
#' an error. In this case, we replace the 'near-zero' eigenvalues of the
#' covariance matrix with non-zero values. In particular, we replace the
#' eigenvalues below the tolerance level, \code{tol}, with the tolerance level.
#' For example, the default value of \code{tol} is 0.001. If we encounter an
#' eigenvalue that is near zero, we will replace it with 0.001 if the covariance
#' matrix estimator is singular.
#'
#' As discussed above, the covariance matrix estimator for each class can be
#' singular, which implies that its determinant is zero. In this case, we also
#' apply the aforementioned eigenvalue adjustment before computing the
#' determinant.
#' 
#' @export
#' @param obj object of type \code{rda} that contains the trained RDA classifier
#' @param newdata matrix containing the unlabeled observations to classify. Each
#' row corresponds to a new observation.
#' @return list with predicted class and discriminant scores for each of the K
#' classes
predict.rda <- function(obj, newdata) {
  if(is.vector(newdata)) {
    newdata <- matrix(newdata, nrow = 1)
  }
  newdata <- as.matrix(newdata)
  
  scores <- sapply(obj$est, function(class_est) {
    cov_inv <- with(class_est$cov_eigen,
                    vectors %*% tcrossprod(diag(1 / values), vectors))
    log_det <- sum(log(class_est$cov_eigen$values))
    apply(newdata, 1, function(obs) {
      x_centered <- as.vector(obs - class_est$xbar)
      obs_quadform <- quadform(A = cov_inv, x = x_centered)
      score <- obs_quadform + log_det - 2 * log(class_est$prior)
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

#' Helper function to optimize the RDA classifier via cross-validation (cv).
#'
#' @export
#' @param x matrix containing the training data. The rows are the sample
#' observations, and the columns are the features.
#' @param y vector of class labels for each training observation
#' @param ... Additional arguments passed to \code{\link{rda}}.
#' @param num_folds the number of cross-validation folds.
#' @param num_lambda The number of values of \code{lambda} to consider
#' @param num_gamma The number of values of \code{gamma} to consider
#' @return list containing the RDA model that minimizes cross-validation as well
#' as a \code{data.frame} that summarizes the cross-validation results.
rda_cv <- function(x, y, ..., num_folds = 10, num_lambda = 5, num_gamma = 5) {
  cv_folds <- cv_partition(y = y, num_folds = num_folds)

  # The grid of tuning parameters
  # We sort the grid by model preference. Models with lambda = 0 are most
  # preferred because they assume that the covariance matrices are equal.
  # This preference is further based on model parsimony.
  tuning_grid <- expand.grid(lambda = seq(0, 1, length = num_lambda),
                             gamma = seq(0, 1, length = num_gamma))
  tuning_grid <- tuning_grid[do.call(order, tuning_grid), ]

  cv_errors <- sapply(cv_folds, function(fold) {
    train_x <- x[fold$training, ]
    train_y <- y[fold$training]
    test_x <- x[fold$test, ]
    test_y <- y[fold$test]

    # TODO: Calculate eigenvalue decomposition
    # TODO: Rewrite 'rda' to accept eigenvalue decomposition
    # TODO: Use this comment for this action once implemented.
    # Rather than recalculating the eigenvalue decomposition for each pair of
    # tuning parameters in the grid, we do this calculation only when the value
    # of lambda changes. This shortcut can reduce the number of computations
    # substantially.
    # if (i %/% num_lambda == 1) { }
    cv_errors <- sapply(seq_len(nrow(tuning_grid)), function(i) {
      rda_out <- with(tuning_grid, rda(x = train_x, y = train_y, ...,
                                       lambda = lambda[i], gamma = gamma[i]))
      sum(predict(rda_out, test_x)$class != test_y)
    })
  })
  cv_errors <- rowSums(cv_errors)

  cv_summary <- cbind(tuning_grid, cv_errors)

  optimal <- which.min(cv_errors)
  lambda <- tuning_grid$lambda[optimal]
  gamma <- tuning_grid$gamma[optimal]
  list(lambda = lambda, gamma = gamma, cv_summary = cv_summary)
}
