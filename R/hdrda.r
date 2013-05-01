#' High-Dimensional Regularized Discriminant Analysis (HDRDA)
#'
#' Given a set of training data, this function builds the HDRDA classifier from
#' Ramey, Stein, and Young (2013). The HDRDA classifier generalizes and improves
#' upon the RDA classifier of Friedman (1989), incorporating a convex combination
#' of the covariance-matrix estimators used in the Linear Discriminant
#' Analysis (LDA) and Quadratic Discriminant Analysis (QDA) classifiers.
#'
#' The matrix of training observations are given in \code{x}. The rows of
#' \code{x} contain the sample observations, and the columns contain the features
#' for each training observation. The vector of class labels given in \code{y}
#' are coerced to a \code{factor}. The length of \code{y} should match the number
#' of rows in \code{x}.
#'
#' The vector, \code{prior}, contains the \emph{a priori} class membership for
#' each class. If \code{prior} is NULL (default), the class membership
#' probabilities are estimated as the sample proportion of observations belonging
#' to each class. Otherwise, \code{prior} should be a vector with the same length
#' as the number of classes in \code{y}. The \code{prior} probabilties should be
#' nonnegative and sum to one. The order of the prior probabilties is assumed to
#' match the levels of factor(y).
#'
#' TODO: Discuss alpha
#' TODO: Discuss lambda
#' TODO: Discuss gamma
#'
#' @export
#' @references Ramey, J. A., Stein, C. K., and Young, D. M. (2013), "A
#' Generalization of Regularized Discriminant Analysis for High-Dimensional
#' Classification."
#' @references Friedman, J. H. (1989), "Regularized Discriminant Analysis,"
#' Journal of American Statistical Association, 84, 405, 165-175.
#' \url{http://www.jstor.org/pss/2289860} (Requires full-text access).
#' @param x matrix containing the training data. The rows are the sample
#' observations, and the columns are the features.
#' @param y vector of class labels for each training observation
#' @param lambda the HDRDA pooling parameter. Must be between 0 and 1,
#' inclusively.
#' @param alpha numeric vector of length \code{K} that scales the convex
#' combination of covariance matrices in the HDRDA classifier. By default, each
#' value is 1. See Ramey et al. (2013) for details.
#' @param gamma a numeric values used for the shrinkage parameter.
#' @param prior vector with prior probabilities for each class. If \code{NULL}
#' (default), then the sample proportion of observations belonging to each class
#' equal probabilities are used. See details.
#' @param tol a threshold for determining nonzero eigenvalues.
#' @return \code{hdrda} object that contains the trained HDRDA classifier
hdrda <- function(x, ...) {
  UseMethod("hdrda")
}

#' @rdname hdrda
#' @method hdrda default
#' @S3method hdrda default
hdrda.default <- function(x, y, lambda = 1, alpha = rep(1, K), gamma = 0,
                         prior = NULL, tol = 1e-6) {
  x <- as.matrix(x)
  y <- as.factor(y)
  lambda <- as.numeric(lambda)
  gamma <- as.numeric(gamma)

  # If p < n, throw warning because it is untested and possibly does not work
  if (ncol(x) < nrow(x)) {
    warning("The sample size exceeds the number of features. This is an untested scenario and may not work currently.")
  }

  if (lambda < 0 || lambda > 1) {
    stop("The value for 'lambda' must be between 0 and 1, inclusively.")
  }

  if (gamma < 0) {
    stop("The value for 'gamma' must be nonnegative.")
  }

  obj <- regdiscrim_estimates(x = x, y = y, cov = FALSE, prior = prior)
  x_centered <- center_data(x = x, y = y)
  K <- obj$num_groups

  if (length(alpha) != K) {
    stop("The length of 'alpha' must equal the number of classes in 'y'.")
  }

  # Computes the eigenvalue decomposition of the pooled sample covariance matrix
  # using the Fast SVD approach.
  cov_pool_eigen <- cov_eigen(x = x, y = y, pool = TRUE, fast = TRUE, tol = tol)

  obj$D_q <- cov_pool_eigen$values
  obj$U_1 <- cov_pool_eigen$vectors
  obj$q <- length(obj$D_q)

  # For each class, we calculate the following quantities necessary to train the
  # HDRDA classifier.
  #   1. \Gamma_k
  #   2. Q_k
  #   3. W_k^{-1}
  for (k in seq_len(K)) {
    # Although 'Gamma_k' is a diagonal matrix, we store only its diagonal
    # elements.
    Gamma <- alpha[k] * lambda * obj$D_q + gamma
    Gamma_inv <- Gamma^(-1)
    X_k <- x_centered[y == levels(y)[k], ]
    n_k <- nrow(X_k)

    XU <- X_k %*% obj$U_1
    Q <- diag(n_k) + alpha[k] * (1 - lambda) * XU %*% tcrossprod(diag(Gamma_inv), XU)

    W_inv <- alpha[k] * (1 - lambda) * diag(Gamma_inv) %*%
      crossprod(XU, solve(Q, XU)) %*% diag(Gamma_inv)
    W_inv <- diag(Gamma_inv) - W_inv

    obj$est[[k]]$n_k <- n_k
    obj$est[[k]]$alpha <- alpha[k]
    obj$est[[k]]$XU <- XU
    obj$est[[k]]$Gamma <- Gamma
    obj$est[[k]]$Q <- Q
    obj$est[[k]]$W_inv <- W_inv
  }

  # Creates an object of type 'hdrda' and adds the 'match.call' to the object
  obj$call <- match.call()
  class(obj) <- "hdrda"
	
  obj
}

#' @param formula A formula of the form \code{groups ~ x1 + x2 + ...} That is,
#' the response is the grouping factor and the right hand side specifies the
#' feature vectors.
#' @param data data frame from which variables specified in \code{formula} are
#' preferentially to be taken.
#' @param ... arguments passed from the \code{formula} to the \code{default}
#' method
#' @rdname hdrda
#' @method hdrda formula
#' @S3method hdrda formula
hdrda.formula <- function(formula, data, ...) {
  # The formula interface includes an intercept. If the user includes the
  # intercept in the model, it should be removed. Otherwise, errors and doom
  # happen.
  # To remove the intercept, we update the formula, like so:
  # (NOTE: The terms must be collected in case the dot (.) notation is used)
  formula <- no_intercept(formula, data)
  
  mf <- model.frame(formula = formula, data = data)
  x <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)

  est <- hdrda.default(x = x, y = y, ...)
  est$call <- match.call()
  est$formula <- formula
  est
}

#' Predicts the class membership of a matrix of unlabeled observations with a
#' trained HDRDA classifier.
#'
#' For a given \code{hdrda} object, we predict the class of each observation (row)
#' of the the matrix given in \code{newdata}.
#'
#' @export
#' @param obj object of type \code{hdrda} that contains the trained HDRDA classifier
#' @param newdata matrix containing the unlabeled observations to classify. Each
#' row corresponds to a new observation.
#' @return list with predicted class and discriminant scores for each of the K
#' classes
predict.hdrda <- function(obj, newdata) {
  if (is.vector(newdata)) {
    newdata <- matrix(newdata, nrow = 1)
  }
  newdata <- as.matrix(newdata)

  scores <- sapply(obj$est, function(class_est) {
    # The call to 'as.vector' removes the attributes returned by 'determinant'
    log_det <- as.vector(determinant(class_est$Q, log = TRUE)$modulus)
    log_det <- log_det + sum(log(class_est$Gamma))

    # Center the 'newdata' by the class sample mean
    x_centered <- scale(newdata, center = class_est$xbar, scale = FALSE)

    # We calculate the quadratic form explicitly for each observation to
    # prevent storing a large 'p x p' inverse matrix in memory. However, note
    # that our approach below increases the number of computations that must be
    # performed for each observation. For the p >> n case, this hardly matters
    # though.
    # The quadratic forms lie on the diagonal of the resulting matrix
    U1_x <- crossprod(obj$U_1, t(x_centered))
    quad_forms <- diag(quadform(class_est$W_inv, U1_x))
    quad_forms + log_det - 2 * log(class_est$prior)
  })

  if (is.vector(scores)) {
    min_scores <- which.min(scores)
    posterior <- exp(-(scores - min(scores)))
  } else {
    min_scores <- apply(scores, 1, which.min)
    # Grabbed code to calculate 'posterior' from MASS:::predict.qda, which
    # handles numerical overflow unlike the more direct:
    # exp(scores) / (1 + exp(sum(scores)))
    posterior <- exp(-(scores - apply(scores, 1, min)))
  }

  class <- with(obj, factor(groups[min_scores], levels = groups))

  list(class = class, scores = scores, posterior = posterior)
}

#' Helper function to optimize the HDRDA classifier via cross-validation (cv).
#'
#' @export
#' @param x matrix containing the training data. The rows are the sample
#' observations, and the columns are the features.
#' @param y vector of class labels for each training observation
#' @param num_folds the number of cross-validation folds.
#' @param num_lambda The number of values of \code{lambda} to consider
#' @param num_gamma The number of values of \code{gamma} to consider
#' @param ... Additional arguments passed to \code{\link{hdrda}}.
#' @return list containing the HDRDA model that minimizes cross-validation as
#' well as a \code{data.frame} that summarizes the cross-validation results.
hdrda_cv <- function(x, y, num_folds = 10, num_lambda = 21, num_gamma = 7, ...) {
  x <- as.matrix(x)
  y <- as.factor(y)
  
  cv_folds <- cv_partition(y = y, num_folds = num_folds)

  # The grid of tuning parameters
  # We sort the grid by model preference. Models with lambda = 0 are most
  # preferred because they assume that the covariance matrices are equal.
  # This preference is further based on model parsimony.
  seq_lambda <- seq(0, 1, length = num_lambda)
  seq_gamma <- 10^seq.int(-1, num_gamma - 2)

  tuning_grid <- expand.grid(lambda = seq_lambda, gamma = seq_gamma)
  tuning_grid <- tuning_grid[do.call(order, tuning_grid), ]

  cv_errors <- sapply(cv_folds, function(fold) {
    train_x <- x[fold$training, ]
    train_y <- y[fold$training]
    test_x <- x[fold$test, ]
    test_y <- y[fold$test]

    hdrda_out <- hdrda(x = train_x, y = train_y, lambda = 1, gamma = 0, ...)

    # For each value of lambda and gamma, we train the HDRDA classifier, classify
    # the test observations, and record the number of classification errors.
    fold_errors <- mapply(function(lambda, gamma) {
      errors <- try({
        # Updates Gamma, Q, and W_inv for each class in hdrda_out
        # If an error is thrown, we return 'NA'.
        hdrda_updated <- update_hdrda(hdrda_out, lambda, gamma)
        sum(predict(hdrda_updated, test_x)$class != test_y)
      }, silent = TRUE)

      errors
    }, tuning_grid$lambda, tuning_grid$gamma)
    fold_errors
  })
  cv_errors <- rowSums(cv_errors)

  cv_summary <- cbind(tuning_grid, cv_errors)

  optimal <- which.min(cv_errors)
  lambda <- tuning_grid$lambda[optimal]
  gamma <- tuning_grid$gamma[optimal]
  list(lambda = lambda, gamma = gamma, cv_summary = cv_summary)
}

#' Helper function to update tuning parameters for the HDRDA classifier
#'
#' TODO
#'
#' @param obj a \code{hdrda} object
#' @param lambda a numeric value between 0 and 1, inclusively
#' @param gamma a numeric value (nonnegative)
#' @return a \code{hdrda} object with updated estimates
update_hdrda <- function(obj, lambda = 1, gamma = 0) {
  for (k in seq_len(obj$num_groups)) {
    XU <- obj$est[[k]]$XU
    alpha <- obj$est[[k]]$alpha
    n_k <- obj$est[[k]]$n_k
    
    Gamma <- alpha * lambda * obj$D_q + gamma
    Gamma_inv <- Gamma^(-1)
    
    Q <- diag(n_k) + alpha * (1 - lambda) * XU %*% tcrossprod(diag(Gamma_inv), XU)
    
    W_inv <- alpha * (1 - lambda) * diag(Gamma_inv) %*%
      crossprod(XU, solve(Q, XU)) %*% diag(Gamma_inv)
    W_inv <- diag(Gamma_inv) - W_inv
    
    obj$est[[k]]$Gamma <- Gamma
    obj$est[[k]]$Q <- Q
    obj$est[[k]]$W_inv <- W_inv
  }
  obj
}
