#' Shrinkage-based Diagonal Quadratic Discriminant Analysis (SDQDA)
#'
#' Given a set of training data, this function builds the DQDA classifier,
#' which is often attributed to Dudoit et al. (2002). To improve the estimation
#' of the pooled variances, we  use a shrinkage method from Pang et al. (2009).
#' This yields the SDQDA classifier.
#'
# The DQDA classifier is a modification to LDA, where the off-diagonal elements
# of the pooled sample covariance matrix are set to zero.
#' 
#' @export
#'
#' @param x training data in matrix form.
#' @param y labels of the training data.
#' @param num_alphas The number of values used to find the optimal amount of shrinkage.
#'
#' @references Dudoit, S., Fridlyand, J., & Speed, T. P. (2002). "Comparison of Discrimination Methods for the Classification of Tumors Using Gene Expression Data," Journal of the American Statistical Association, 97, 457, 77-87. 
#' @references Pang, H., Tong, T., & Zhao, H. (2009). "Shrinkage-based Diagonal Discriminant Analysis and Its Applications in High-Dimensional Data," Biometrics, 65, 4, 1021-1029.
#' @return sdqda obj
sdqda <- function(x, y, num_alphas = 101) {
	obj <- list()
	y <- factor(y)
	obj$labels <- y
	obj$N <- length(y)
	obj$p <- ncol(x)
	obj$groups <- levels(y)
	obj$num_groups <- nlevels(obj$labels)
	
	obj$est <- foreach(k=levels(y)) %do% {
		stats <- list()
		x_k <- x[which(y == k), ]
		stats$label <- k
		stats$n_k <- n_k <- nrow(x_k)
		stats$xbar <- colMeans(x_k)
		stats$var <- (n_k - 1) * apply(x_k, 2, var) / n_k
		stats$var_shrink <- var_shrinkage(
			N = n_k,
			K = 1,
			var_feature = stats$var,
			num_alphas = num_alphas,
			t = -1
		)
		stats
	}

	class(obj) <- "sdqda"
	obj	
}

#' SDQDA prediction of the class membership of a matrix of new observations.
#'
# The DQDA classifier is a modification to LDA, where the off-diagonal elements
# of the pooled sample covariance matrix are set to zero. To improve the estimation
#' of the pooled variances, we use a shrinkage method from Pang et al. (2009).
#' This yields the SDQDA classifier.
#' 
#' @export
#'
#' @param obj trained SDQDA object
#' @param newdata matrix of observations to predict. Each row corresponds to a new observation.
#'
#' @references Dudoit, S., Fridlyand, J., & Speed, T. P. (2002). "Comparison of Discrimination Methods for the Classification of Tumors Using Gene Expression Data," Journal of the American Statistical Association, 97, 457, 77-87. 
#' @references Pang, H., Tong, T., & Zhao, H. (2009). "Shrinkage-based Diagonal Discriminant Analysis and Its Applications in High-Dimensional Data," Biometrics, 65, 4, 1021-1029.
#' @return list predicted class memberships of each row in newdata
predict_sdqda <- function(obj, newdata) {
	if (!inherits(obj, "sdqda"))  {
		stop("obj not of class 'sdqda'")
	}
	if(is.vector(newdata)) newdata <- matrix(newdata, nrow = 1)

	scores <- apply(newdata, 1, function(obs) {
		sapply(obj$est, function(class_est) {
			sum((obs - class_est$xbar)^2 * class_est$var_shrink) 
		})
	})
	
	if(is.vector(scores)) {
		min_scores <- which.min(scores)
	} else {
		min_scores <- apply(scores, 2, which.min)
	}

	predicted <- factor(obj$groups[min_scores], levels = obj$groups)
	
	list(scores = scores, predicted = predicted)
}
