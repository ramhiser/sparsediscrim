#' Computes estimates and ancillary information for diagonal classifiers
#'
#' Computes the maximum likelihood estimators (MLEs) for each class under the
#' assumption of multivariate normality for each class. Also, computes ancillary
#' information necessary for classifier summary, such as sample size, the number
#' of features, etc.
#'
#' This function computes the common estimates and ancillary information used in
#' all of the diagonal classifiers in the \code{diagdiscrim} package.
#'
#' The matrix of training observations are given in \code{x}. The rows of \code{x}
#' contain the sample observations, and the columns contain the features for each
#' training observation.
#'
#' The vector of class labels given in \code{y} are coerced to a \code{factor}.
#' The length of \code{y} should match the number of rows in \code{x}.
#'
#' An error is thrown if a given class has less than 2 observations because the
#' variance for each feature within a class cannot be estimated with less than 2
#' observations.
#'
#' The vector, \code{prior}, contains the \emph{a priori} class membership for
#' each class. If \code{prior} is NULL (default), the class membership
#' probabilities are estimated as the sample proportion of observations belonging
#' to each class. Otherwise, \code{prior} should be a vector with the same length
#' as the number of classes in \code{y}. The \code{prior} probabilties should be
#' nonnegative and sum to one.
#' 
#' @param matrix containing the training data. The rows are the sample
#' observations, and the columns are the features.
#' @param y vector of class labels for each training observation
#' @param prior vector with prior probabilities for each class. If NULL
#' (default), then equal probabilities are used. See details.
#' @param pool logical value. If TRUE, calculates the pooled sample variances
#' for each class.
#' @return list of MLEs and necessary ancillary information
diag_estimates <- function(x, y, prior = NULL, pool = FALSE) {
  obj <- list()
	obj$labels <- y
	obj$N <- length(y)
	obj$p <- ncol(x)
	obj$groups <- levels(y)
	obj$num_groups <- nlevels(y)

  # Error Checking
  if (!is.null(prior)) {
    if (length(prior) != obj$num_groups) {
      stop("The number of 'prior' probabilities must match the number of classes in 'y'.")
    }
    if (any(prior <= 0)) {
      stop("The 'prior' probabilities must be nonnegative.")
    }
    if (sum(prior) != 1) {
      stop("The 'prior' probabilities must sum to one.")
    }
  }
  if (any(table(y) < 2)) {
    stop("There must be at least 2 observations in each class.")
  }

  # By default, we estimate the 'a priori' probabilties of class membership with
  # the MLEs (the sample proportions).
  if (is.null(prior)) {
    prior <- as.vector(table(y) / length(y))
  }

  # For each class, we calculate the MLEs for each parameter used in the DLDA
  # classifier. The 'est' list contains the estimators for each class.
  obj$est <- tapply(seq_along(y), y, function(i) {
    stats <- list()
    stats$n <- length(i)
    stats$xbar <- colMeans(x[i,])
    stats$var <- with(stats, (n - 1) * apply(x[i,], 2, var) / n)
    stats
  })

  # Calculates the pooled variance across all classes.
  if (pool) {
    obj$var_pool <- Reduce('+', lapply(obj$est, function(x) x$n * x$var)) / obj$N
    obj$est <- lapply(obj$est, function(class_est) {
      class_est$var <- obj$var_pool
      class_est
    })
  }

  # Add each element in 'prior' to the corresponding obj$est$prior
  for(k in seq_len(obj$num_groups)) {
    obj$est[[k]]$prior <- prior[k]
  }
  obj
}
