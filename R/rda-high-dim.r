#' High-Dimensional Regularized Discriminant Analysis (HDRDA)
#'
#' Given a set of training data, this function builds the HDRDA classifier from
#' Ramey, Stein, and Young (2017). Specially designed for small-sample,
#' high-dimensional data, the HDRDA classifier incorporates dimension reduction
#' and covariance-matrix shrinkage to enable a computationally efficient
#' classifier.
#'
#' The HDRDA classifier utilizes a covariance-matrix estimator that is a convex
#' combination of the covariance-matrix estimators used in the Linear
#' Discriminant Analysis (LDA) and Quadratic Discriminant Analysis (QDA)
#' classifiers. For each of the \code{K} classes given in \code{y},
#' \eqn{(k = 1, \ldots, K)}, we first define this convex combination as
#' \deqn{\hat{\Sigma}_k(\lambda) = (1 - \lambda) \hat{\Sigma}_k
#' + \lambda \hat{\Sigma},}
#' where \eqn{\lambda \in [0, 1]} is the \emph{pooling} parameter. We then
#' calculate the covariance-matrix estimator
#' \deqn{\tilde{\Sigma}_k = \alpha_k \hat{\Sigma}_k(\lambda) + \gamma I_p,}
#' where \eqn{I_p} is the \eqn{p \times p} identity matrix. The matrix
#' \eqn{\tilde{\Sigma}_k} is substituted into the HDRDA classifier. See Ramey et
#' al. (2017) for more details.
#'
#' The matrix of training observations are given in \code{x}. The rows of
#' \code{x} contain the sample observations, and the columns contain the features
#' for each training observation. The vector of class labels given in \code{y}
#' are coerced to a \code{factor}. The length of \code{y} should match the number
#' of rows in \code{x}.
#'
#' The vector \code{prior} contains the \emph{a priori} class membership for
#' each class. If \code{prior} is \code{NULL} (default), the class membership
#' probabilities are estimated as the sample proportion of observations
#' belonging to each class. Otherwise, \code{prior} should be a vector with the
#' same length as the number of classes in \code{y}. The \code{prior}
#' probabilities should be nonnegative and sum to one. The order of the prior
#' probabilities is assumed to match the levels of \code{factor(y)}.
#'
#' @export
#' @references Ramey, J. A., Stein, C. K., and Young, D. M. (2017),
#' "High-Dimensional Regularized Discriminant Analysis."
#' \url{https://arxiv.org/abs/1602.01182}.
#' @references Friedman, J. H. (1989), "Regularized Discriminant Analysis,"
#' Journal of American Statistical Association, 84, 405, 165-175.
#' \url{http://www.jstor.org/pss/2289860} (Requires full-text access).
#' @param x matrix containing the training data. The rows are the sample
#' observations, and the columns are the features.
#' @param y vector of class labels for each training observation
#' @param lambda the HDRDA pooling parameter. Must be between 0 and 1,
#' inclusively.
#' @param gamma a numeric values used for the shrinkage parameter.
#' @param shrinkage_type the type of covariance-matrix shrinkage to apply. By
#' default, a ridge-like shrinkage is applied. If \code{convex} is given, then
#' shrinkage similar to Friedman (1989) is applied. See Ramey et al. (2017) for
#' details.
#' @param prior vector with prior probabilities for each class. If \code{NULL}
#' (default), then the sample proportion of observations belonging to each class
#' equal probabilities are used. See details.
#' @param tol a threshold for determining nonzero eigenvalues.
#' @return \code{rda_high_dim} object that contains the trained HDRDA classifier
rda_high_dim <- function(x, ...) {
  UseMethod("rda_high_dim")
}

#' @rdname rda_high_dim
#' @export
rda_high_dim.default <- function(x, y, lambda = 1, gamma = 0,
                                 shrinkage_type = c("ridge", "convex"), prior = NULL,
                                 tol = 1e-6, ...) {
  x <- as.matrix(x)
  y <- as.factor(y)
  lambda <- as.numeric(lambda)
  gamma <- as.numeric(gamma)
  shrinkage_type <- match.arg(shrinkage_type)
  
  if (lambda < 0 || lambda > 1) {
    stop("The value for 'lambda' must be between 0 and 1, inclusively.")
  }
  
  if (gamma < 0) {
    stop("The value for 'gamma' must be nonnegative.")
  }
  
  obj <- regdiscrim_estimates(x = x, y = y, cov = FALSE, prior = prior, ...)
  x_centered <- center_data(x = x, y = y)
  K <- obj$num_groups
  
  # Computes the eigenvalue decomposition of the pooled sample covariance matrix
  # using the Fast SVD approach.
  cov_pool_eigen <- cov_eigen(x = x, y = y, pool = TRUE, fast = TRUE, tol = tol)
  
  obj$D_q <- cov_pool_eigen$values
  obj$U1 <- cov_pool_eigen$vectors
  obj$q <- length(obj$D_q)
  obj$shrinkage_type <- shrinkage_type
  obj$lambda <- lambda
  obj$gamma <- gamma
  
  if (shrinkage_type == "ridge") {
    alpha <- 1
  } else {
    # shinkage_family == "convex"
    alpha <- 1 - gamma
  }
  
  # Transforms the centered data
  XU <- x_centered %*% obj$U1
  
  # For each class, we calculate the following quantities necessary to train the
  # HDRDA classifier.
  #   1. \Gamma_k
  #   2. Q_k
  #   3. W_k^{-1}
  for (k in seq_len(K)) {
    X_k <- x_centered[y == levels(y)[k], , drop=FALSE]
    n_k <- nrow(X_k)
    
    # Extracts the transformed, centered data
    # No need to calculate it for the classes individually
    XU_k <- XU[y == levels(y)[k], , drop=FALSE]
    
    # Transforms the sample mean to the lower dimension
    xbar_U1 <- crossprod(obj$U1, obj$est[[k]]$xbar)
    
    # In the special case that (lambda, gamma) = (0, 0), the improvements to
    # HDRDA's speed via the Sherman-Woodbury formula are not applicable because
    # Gamma = 0. In this case, we calculate W_inv directly. If the matrix is
    # singular, a slight amount of shrinkage is applied.
    if (lambda == 0 && gamma == 0) {
      W_k <- cov_mle(XU_k)
      W_inv <- try(solve_chol(W_k), silent=TRUE)
      
      if (inherits(W_inv, "try-error")) {
        W_k <- W_k + diag(0.001, nrow=nrow(W_k), ncol=ncol(W_k))
        W_inv <- solve_chol(W_k)
      }
      
      Gamma <- matrix(0, nrow=obj$q, ncol=obj$q)
      Q <- diag(n_k)
    } else {
      # Although 'Gamma_k' is a diagonal matrix, we store only its diagonal
      # elements.
      Gamma <- alpha * lambda * obj$D_q + gamma
      Gamma_inv <- diag(Gamma^(-1))
      
      # X_k %*% U_1 %*% Gamma^{-1} is computed repeatedly in the equations.
      # We compute the matrix once and use it where necessary to avoid
      # unnecessary computations.
      XU_Gamma_inv <- XU_k %*% Gamma_inv
      Q <- diag(n_k) + alpha * (1 - lambda) / n_k * tcrossprod(XU_Gamma_inv, XU_k)
      W_inv <- alpha * (1 - lambda) / n_k * crossprod(XU_Gamma_inv, solve(Q, XU_Gamma_inv))
      W_inv <- Gamma_inv - W_inv
    }
    
    obj$est[[k]]$n_k <- n_k
    obj$est[[k]]$alpha <- alpha
    obj$est[[k]]$XU <- XU_k
    obj$est[[k]]$xbar_U1 <- xbar_U1
    obj$est[[k]]$Gamma <- Gamma
    obj$est[[k]]$Q <- Q
    obj$est[[k]]$W_inv <- W_inv
  }
  
  # Creates an object of type 'rda_high_dim' and adds the 'match.call' to the object
  obj$call <- match.call()
  class(obj) <- "rda_high_dim"
  
  obj
}

#' @param formula A formula of the form \code{groups ~ x1 + x2 + ...} That is,
#' the response is the grouping factor and the right hand side specifies the
#' feature vectors.
#' @param data data frame from which variables specified in \code{formula} are
#' preferentially to be taken.
#' @param ... arguments passed from the \code{formula} to the \code{default}
#' method
#' @rdname rda_high_dim
#' @importFrom stats model.frame model.matrix model.response
#' @export
rda_high_dim.formula <- function(formula, data, ...) {
  # The formula interface includes an intercept. If the user includes the
  # intercept in the model, it should be removed. Otherwise, errors and doom
  # happen.
  # To remove the intercept, we update the formula, like so:
  # (NOTE: The terms must be collected in case the dot (.) notation is used)
  formula <- no_intercept(formula, data)
  
  mf <- model.frame(formula = formula, data = data)
  x <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)
  
  est <- rda_high_dim.default(x = x, y = y, ...)
  est$call <- match.call()
  est$formula <- formula
  est
}

#' Outputs the summary for a HDRDA classifier object.
#'
#' Summarizes the trained rda_high_dim classifier in a nice manner.
#'
#' @param x object to print
#' @param ... unused
#' @export
print.rda_high_dim <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("Sample Size:\n")
  print(x$N)
  cat("Number of Features:\n")
  print(x$p)
  cat("Classes:\n")
  print(x$groups)
  cat("Prior Probabilities:\n")
  print(sapply(x$est, function(z) z$prior))
  cat("Reduced dimension, q:\n")
  print(x$q)
  cat("Shrinkage type:\n")
  print(x$shrinkage_type)
  cat("Lambda:\n")
  print(x$lambda)
  cat("Gamma:\n")
  print(x$gamma)
}

#' Predicts the class membership of a matrix of unlabeled observations with a
#' trained HDRDA classifier.
#'
#' For a given \code{rda_high_dim} object, we predict the class of each observation
#' (row) of the the matrix given in \code{newdata}.
#'
#' @rdname rda_high_dim
#' @export
#' @param object object of type \code{rda_high_dim} that contains the trained HDRDA
#' classifier
#' @param newdata matrix containing the unlabeled observations to classify. Each
#' row corresponds to a new observation.
#' @param projected logical indicating whether \code{newdata} have already been
#' projected to a q-dimensional subspace. This argument can yield large gains in
#' speed when the linear transformation has already been performed.
#' @return list with predicted class and discriminant scores for each of the K
#' classes
predict.rda_high_dim <- function(object, newdata, projected = FALSE, ...) {
  newdata <- as.matrix(newdata)
  
  scores <- sapply(object$est, function(class_est) {
    if (object$lambda == 0 && object$gamma == 0) {
      # Want: log(det(W_k)) = -log(det(W_k_inv))
      log_det <- -log_determinant(class_est$W_inv)
    } else {
      log_det <- log_determinant(class_est$Q)
    }
    
    if (projected) {
      # The newdata have already been projected. Yay for speed!
      # Center the 'newdata' by the projected class sample mean
      U1_x <- scale(newdata, center = class_est$xbar_U1, scale = FALSE)
      
      quad_forms <- diag(drop(tcrossprod(U1_x %*% class_est$W_inv, U1_x)))
    } else {
      # Center the 'newdata' by the class sample mean
      x_centered <- scale(newdata, center = class_est$xbar, scale = FALSE)
      
      U1_x <- crossprod(object$U1, t(x_centered))
      
      quad_forms <- apply(U1_x, 2, function(z) {
        quadform(class_est$W_inv, z)
      })
    }
    
    quad_forms + log_det - 2 * log(class_est$prior)
  })
  
  if (is.vector(scores)) {
    # When sapply above returns a vector, the naming is thrown off.
    # Hence, we rename it to the groups.
    names(scores) <- object$groups
    
    min_scores <- which.min(scores)
    posterior <- exp(-(scores - min(scores)))
    posterior <- posterior / sum(posterior)
  } else {
    min_scores <- apply(scores, 1, which.min)
    # Grabbed code to calculate 'posterior' from MASS:::predict.qda, which
    # handles numerical overflow unlike the more direct:
    # exp(scores) / (1 + exp(sum(scores)))
    posterior <- exp(-(scores - apply(scores, 1, min)))
    posterior <- posterior / rowSums(posterior)
  }
  
  class <- with(object, factor(groups[min_scores], levels = groups))
  
  list(class = class, scores = scores, posterior = posterior)
}

#' Helper function to optimize the HDRDA classifier via cross-validation
#'
#' For a given data set, we apply cross-validation (cv) to select the optimal
#' HDRDA tuning parameters.
#'
#' The number of cross-validation folds is given in \code{num_folds}.
#'
#' @export
#' @importFrom dplyr arrange
#' @param x matrix containing the training data. The rows are the sample
#' observations, and the columns are the features.
#' @param y vector of class labels for each training observation
#' @param num_folds the number of cross-validation folds.
#' @param num_lambda The number of values of \code{lambda} to consider
#' @param num_gamma The number of values of \code{gamma} to consider
#' @param shrinkage_type the type of covariance-matrix shrinkage to apply. By
#' default, a ridge-like shrinkage is applied. If \code{convex} is given, then
#' shrinkage similar to Friedman (1989) is applied. See Ramey et al. (2017) for
#' details.
#' @param verbose If set to \code{TRUE}, summary information will be outputted
#' as the optimal model is being determined.
#' @param ... Additional arguments passed to \code{\link{rda_high_dim}}.
#' @return list containing the HDRDA model that minimizes cross-validation as
#' well as a \code{data.frame} that summarizes the cross-validation results.
rda_high_dim_cv <- function(x, y, num_folds = 10, num_lambda = 21, num_gamma = 8,
                            shrinkage_type=c("ridge", "convex"), verbose=FALSE, ...) {
  x <- as.matrix(x)
  y <- as.factor(y)
  shrinkage_type <- match.arg(shrinkage_type)
  
  cv_folds <- cv_partition(y = y, num_folds = num_folds)
  
  # The grid of tuning parameters
  # We sort the grid by model preference. Models with lambda = 0 are most
  # preferred because they assume that the covariance matrices are equal. This
  # preference is further based on model parsimony.
  seq_lambda <- seq(0, 1, length = num_lambda)
  
  if (shrinkage_type == "ridge") {
    seq_gamma <- c(0, 10^seq.int(-2, num_gamma - 4))
  } else  {
    # shrinkage_type == "convex"
    seq_gamma <- seq(0, 1, length = num_gamma)
  }
  
  tuning_grid <- expand.grid(lambda = seq_lambda, gamma = seq_gamma)
  tuning_grid <- dplyr::arrange(tuning_grid, lambda, gamma)
  
  cv_errors <- sapply(seq_along(cv_folds), function(fold_i) {
    if (verbose) {
      cat("CV Fold: ", fold_i, " of ", num_folds, "\n")
    }
    fold <- cv_folds[[fold_i]]
    train_x <- x[fold$training, ]
    train_y <- y[fold$training]
    test_x <- x[fold$test, ]
    test_y <- y[fold$test]
    
    hdrda_out <- rda_high_dim(x = train_x, y = train_y, lambda = 1, gamma = 0, ...)
    
    # Projects the test data to the q-dimensional subspace.
    # No need to do this for each lambda/gamma pair.
    # NOTE: It's not centered yet.
    test_x <- test_x %*% hdrda_out$U1
    
    # For each value of lambda and gamma, we train the HDRDA classifier,
    # classify the test observations, and record the number of classification
    # errors.
    fold_errors <- mapply(function(lambda, gamma) {
      errors <- try({
        # Updates Gamma, Q, and W_inv for each class in hdrda_out
        # If an error is thrown, we return 'NA'.
        hdrda_updated <- update_rda_high_dim(hdrda_out, lambda, gamma)
        # Call to predict.rda_high_dim is explicit to pass R CMD CHECK.
        sum(predict.rda_high_dim(hdrda_updated, test_x, projected = TRUE)$class != test_y)
      }, silent = TRUE)
      
      errors
    }, tuning_grid$lambda, tuning_grid$gamma)
    fold_errors
  })
  cv_errors <- rowSums(cv_errors)
  error_rate <- cv_errors / nrow(x)
  
  cv_summary <- cbind(tuning_grid, cv_errors, error_rate)
  
  optimal <- which.min(cv_errors)
  lambda <- tuning_grid$lambda[optimal]
  gamma <- tuning_grid$gamma[optimal]
  
  # Trains a classifier based on optimal model
  hdrda_out <- rda_high_dim(x=x, y=y, lambda=lambda, gamma=gamma,
                            shrinkage_type=shrinkage_type)
  
  # Adds optimal parameters and cv_summary to classifier object
  hdrda_out$lambda <- lambda
  hdrda_out$gamma <- gamma
  hdrda_out$cv_summary <- cv_summary
  class(hdrda_out) <- c("hdrda_cv", "rda_high_dim")
  
  hdrda_out
}

#' Plots a heatmap of cross-validation error grid for a HDRDA classifier object.
#'
#' Uses \code{\link[ggplot2]{ggplot}} to plot a heatmap of the training error
#' grid.
#'
#' @param x object to plot
#' @param ... unused
#' @export
#' @importFrom ggplot2 ggplot aes scale_fill_gradient labs scale_x_discrete
#' @importFrom ggplot2 scale_y_discrete theme labs geom_tile element_blank
plot.rda_high_dim_cv <- function(x, ...) {
  cv_summary <- x$cv_summary
  cv_summary <- within(cv_summary, {
    lambda <- round(lambda, 3)
    gamma <- round(gamma, 3)
  })
  
  # Fixes NOTE from R CMD CHECK
  # "no visible binding for global variable 'error_rate'"
  error_rate <- 1
  
  p <- ggplot(cv_summary, aes(factor(gamma), factor(lambda)))
  p <- p + geom_tile(aes(fill=error_rate), colour="white")
  p <- p + scale_fill_gradient(low="white",
                               high="steelblue",
                               name="CV Error Rate")
  p <- p + labs(x=expression(gamma), y=expression(lambda))
  p <- p + scale_x_discrete(expand=c(0, 0))
  p <- p + scale_y_discrete(expand=c(0, 0))
  p <- p + theme(axis.ticks=element_blank())
  p
}

#' Helper function to update tuning parameters for the HDRDA classifier
#'
#' This function updates some of the quantities in the HDRDA classifier based on
#' updated values of \code{lambda} and \code{gamma}. The update can greatly
#' expedite cross-validation to examine a large grid of values for \code{lambda}
#' and \code{gamma}.
#'
#' @param obj a \code{rda_high_dim} object
#' @param lambda a numeric value between 0 and 1, inclusively
#' @param gamma a numeric value (nonnegative)
#' @return a \code{rda_high_dim} object with updated estimates
update_rda_high_dim <- function(obj, lambda = 1, gamma = 0) {
  # In the special case that (lambda, gamma) = (0, 0), the improvements to
  # HDRDA's speed via the Sherman-Woodbury formula are not applicable because
  # Gamma = 0. In this case, we calculate W_inv directly.
  if (lambda == 0 && gamma == 0) {
    Gamma <- matrix(0, nrow=obj$q, ncol=obj$q)
  } else {
    # NOTE: alpha_k is constant across all classes, so that alpha_k = alpha_1
    # for all k. As a result, Gamma and Gamma_inv are constant across all k. We
    # compute both before looping through all K classes.
    Gamma <- obj$est[[1]]$alpha * lambda * obj$D_q + gamma
    Gamma_inv <- diag(Gamma^(-1))
  }
  
  for (k in seq_len(obj$num_groups)) {
    n_k <- obj$est[[k]]$n_k
    if (lambda == 0 && gamma == 0) {
      Q <- diag(n_k)
      
      W_k <- cov_mle(obj$est[[k]]$XU)
      W_inv <- try(solve_chol(W_k), silent=TRUE)
      if (inherits(W_inv, "try-error")) {
        W_k <- W_k + diag(0.001, nrow=nrow(W_k), ncol=ncol(W_k))
        W_inv <- solve_chol(W_k)
      }
    }
    else {
      # X_k %*% U_1 %*% Gamma^{-1} is computed repeatedly in the equations.
      # We compute the matrix once and use it where necessary to avoid
      # unnecessary computations.
      XU_Gamma_inv <- obj$est[[k]]$XU %*% Gamma_inv
      Q <- diag(n_k) +
        obj$est[[k]]$alpha * (1 - lambda) / n_k * tcrossprod(XU_Gamma_inv, obj$est[[k]]$XU)
      
      W_inv <- obj$est[[k]]$alpha * (1 - lambda) / n_k *
        crossprod(XU_Gamma_inv, solve(Q, XU_Gamma_inv))
      W_inv <- Gamma_inv - W_inv
    }
    obj$est[[k]]$Gamma <- Gamma
    obj$est[[k]]$Q <- Q
    obj$est[[k]]$W_inv <- W_inv
  }
  obj$lambda <- lambda
  obj$gamma <- gamma
  obj
}
