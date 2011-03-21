# Regularized Shrinkage-based Diagonal Discriminant Analysis (RSDDA)
# The RSDDA classifier is a regularized version of SDQDA, where the off-diagonal elements
# of the pooled sample covariance matrix are set to zero.
# We use a shrinkage method based on Pang, Tong, and Zhao (2009).
# We assume the first column is named "labels" and holds a factor vector,
# which contains the class labels.

# We assume the first column is named "labels" and holds a factor vector,
# which contains the class labels.

# num_alphas: the grid size of alphas considered in estimating each shrinkage parameter alpha
rsdda <- function(train_df, num_alphas = 101) {
	obj <- list()
	
	N <- nrow(train_df)
	obj$training <- train_df
	obj$N <- N
	obj$classes <- levels(train_df$labels)
	obj$num_classes <- nlevels(train_df$labels)

	estimators <- dlply(obj$training, .(labels), function(df_k) {
		n_k <- nrow(df_k)
		pi_k <- n_k / N
		
		xbar <- as.vector(colMeans(df_k))
		
		sum_squares <- (n_k - 1) * apply(class.x, 2, function(col) { var(col) })
		var <- sum_squares / n_k
		var_shrink <- var_shrinkage(N = n_k, K = 1, var_feature = var, num_alphas = num_alphas, t = -1)
		
		list(xbar = xbar, var_k = var_shrink, sum_squares = sum_squares, n = n_k, pi_k = pi_k)
	})
	
	var_pool <- colSums(laply(estimators, function(class_est) class_est$sum_squares)) / N
	var_pool_shrink <- var_shrinkage(N = N, K = num.classes, var_feature = var_pool, num_alphas = num_alphas, t = -1)
	
	obj$estimators <- llply(estimators, function(class_estimators) {
		class_estimators$var_pool <- var_pool_shrink
		class_estimators
	})
	
	class(obj) <- "rsdda"

	obj
}

# Estimates the optimal value of lambda used in RSDDA by leave-k-out crossvalidation.
# grid_size: the grid size of lambda_grid considered in estimating the regularization parameter lambda
# k: Leave-k-out crossvalidation is used to estimate the optimal value of lambda.
model.select.rsdda <- function(object, grid_size = 5, k = 1) {
	if (!inherits(object, "rsdda"))  {
		stop("object not of class 'rsdda'")
	}
	# Estimate lambda with Leave-k-out crossvalidation.
	if(k > 1) {
		k <- 1
		warning("Leave-k-out crossvalidation for k > 1 is not implemented yet. Using k = 1 instead...\n")
	}
	
	lambda_grid <- seq(0, 1, length = grid_size)
	
	# Here, we calculate the Leave-K-out Error Rates for each lambda in the grid.
	# NOTE: We are only simultaneously diagonalizing the covariance matrices ONCE.
	# If we did this for each observation left out using crossvalidation, the model
	# selection process would be take too long to run, especially in simulations
	# with a large number of replications.
	predictions <- laply(seq_len(object$N), function(i) {
		loo.obj <- rsdda(object$training[-i, ], num_alphas = grid_size)
		laply(lambda_grid, function(lambda) {
			predict.rsdda(loo.obj, object$training[i, -1], lambda = lambda)
		})
	})

	error.rates.loo <- apply(predictions, 2, function(predictions.lambda) {
		mean(object$training[, 1] != predictions.lambda)
	})
	
	# Find the lambdas that attain the minimum Leave-K-out error rate.
	optimal_lambdas <- lambda_grid[which(error.rates.loo == min(error.rates.loo))]
	
	# Finally, we add the optimal lambda to the RSDDA object.
	# If there is a tie for optimal lambda, randomly choose one from the optimal ones.
	object$lambda <- sample(optimal_lambdas, 1)
	
	return(object)
}

predict.rsdda <- function(object, newdata, num_lambdas = 5, lambda = NULL, verbose = FALSE) {
	if (!inherits(object, "rsdda"))  {
		stop("object not of class 'rsdda'")
	}
	
	if(is.null(lambda) || !is.numeric(lambda)) {
		if(verbose) cat("Model selection for RSDDA\n")
		object <- model.select.rsdda(object, grid_size = num_lambdas)
		if(verbose) cat("Model selection for RSDDA...done!\n")
		if(verbose) cat("RSDDA Lambda Parameter:", object$lambda, "\n")
	}
	else {
		object$lambda <- lambda
	}
	
	newdata <- data.matrix(newdata)

	predictions <- apply(newdata, 1, function(obs) {
		scores <- sapply(object$estimators, function(class_est) {
			var_rsdda <- (class_est$var_k)^(1 - object$lambda) * (class_est$var_pool)^(object$lambda)
			sum((obs - class_est$xbar)^2 * var_rsdda) - sum(log(var_rsdda)) - 2 * log(class_est$pi_k)
		})
		
		prediction <- object$classes[which.min(scores)]
		prediction
	})
	
	predictions <- factor(predictions, levels = object$classes)

	predictions
}