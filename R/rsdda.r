# Regularized Shrinkage-based Diagonal Discriminant Analysis (RSDDA)
# The RSDDA classifier is a regularized version of SDQDA, where the off-diagonal elements
# of the pooled sample covariance matrix are set to zero.
# We use a shrinkage method based on Pang, Tong, and Zhao (2009).
# We assume the first column is named "labels" and holds a factor vector,
# which contains the class labels.

# We assume the first column is named "labels" and holds a factor vector,
# which contains the class labels.

# num_alphas: the grid size of alphas considered in estimating each shrinkage parameter alpha
rsdda <- function(train_df, num_alphas = 5, jointdiag = "none", verbose = FALSE, ...) {
	rsdda.obj <- list()
	rsdda.obj$training <- train_df
	
	if(jointdiag != "none") {
		if(verbose) message("Simultaneously diagonalizing covariance matrices... ", appendLF = FALSE)
		joint.diag.out <- joint.diagonalization(rsdda.obj$training, method = jointdiag)
		rsdda.obj$training <- joint.diag.out$transformed.df
		rsdda.obj$jointdiag.B <- joint.diag.out$B
		rsdda.obj$jointdiag.method <- joint.diag.out$method
		if(verbose) message("done!")
	}

	if(verbose) message("Building RSDDA classifier... ", appendLF = FALSE)
	N <- nrow(rsdda.obj$training)
	num.classes <- nlevels(rsdda.obj$training$labels)
	
	train_x <- as.matrix(rsdda.obj$training[,-1])
	dimnames(train_x) <- NULL
	
	estimators <- dlply(rsdda.obj$training, .(labels), function(df_k) {
		class.x <- as.matrix(df_k[, -1])
		dimnames(class.x) <- NULL
		
		n_k <- nrow(class.x)
		pi_k <- n_k / N
		
		xbar <- as.vector(colMeans(class.x))
		
		sum_squares <- apply(class.x, 2, function(col) {
			(n_k - 1) * var(col)
		})
		
		var <- apply(class.x, 2, function(col) {
			(n_k - 1) * var(col) / n_k
		})
		
		var_shrink <- var_shrinkage(N = n_k, K = 1, var_feature = var, num_alphas = num_alphas, t = -1)
		
		list(xbar = xbar, var.k = var_shrink, sum_squares = sum_squares, n = n_k, pi_k = pi_k)
	})
	
	var_pool <- colSums(laply(estimators, function(class_est) class_est$sum_squares)) / N
	var.pool.shrink <- var_shrinkage(N = N, K = num.classes, var_feature = var_pool, num_alphas = num_alphas, t = -1)
	
	estimators <- llply(estimators, function(class_estimators) {
		class_estimators$var.pool <- var.pool.shrink
		class_estimators
	})
	
	if(verbose) message("done!")
	
	rsdda.obj$N <- N
	rsdda.obj$classes <- levels(rsdda.obj$training$labels)
	rsdda.obj$estimators <- estimators
	
	class(rsdda.obj) <- "rsdda"

	rsdda.obj
}

# Estimates the optimal value of lambda used in RSDDA by leave-k-out crossvalidation.
# grid.size: the grid size of lambda.grid considered in estimating the regularization parameter lambda
# k: Leave-k-out crossvalidation is used to estimate the optimal value of lambda.
model.select.rsdda <- function(object, grid.size = 5, k = 1) {
	if (!inherits(object, "rsdda"))  {
		stop("object not of class 'rsdda'")
	}
	# Estimate lambda with Leave-k-out crossvalidation.
	if(k > 1) {
		k <- 1
		warning("Leave-k-out crossvalidation for k > 1 is not implemented yet. Using k = 1 instead...\n")
	}
	
	lambda.grid <- seq(0, 1, length = grid.size)
	
	# Here, we calculate the Leave-K-out Error Rates for each lambda in the grid.
	# NOTE: We are only simultaneously diagonalizing the covariance matrices ONCE.
	# If we did this for each observation left out using crossvalidation, the model
	# selection process would be take too long to run, especially in simulations
	# with a large number of replications.
	predictions <- laply(seq_len(object$N), function(i) {
		loo.rsdda.obj <- rsdda(object$training[-i, ], num_alphas = grid.size)
		laply(lambda.grid, function(lambda) {
			predict.rsdda(loo.rsdda.obj, object$training[i, -1], lambda = lambda)
		})
	})

	error.rates.loo <- apply(predictions, 2, function(predictions.lambda) {
		mean(object$training[, 1] != predictions.lambda)
	})
	
	# Find the lambdas that attain the minimum Leave-K-out error rate.
	optimal.lambdas <- lambda.grid[which(error.rates.loo == min(error.rates.loo))]
	
	# Finally, we add the optimal lambda to the RSDDA object.
	# If there is a tie for optimal lambda, randomly choose one from the optimal ones.
	object$lambda <- sample(optimal.lambdas, 1)
	
	return(object)
}

predict.rsdda <- function(object, newdata, num.lambdas = 5, lambda = NULL, verbose = FALSE) {
	if (!inherits(object, "rsdda"))  {
		stop("object not of class 'rsdda'")
	}
	
	if(is.null(lambda) || !is.numeric(lambda)) {
		if(verbose) cat("Model selection for RSDDA\n")
		object <- model.select.rsdda(object, grid.size = num.lambdas)
		if(verbose) cat("Model selection for RSDDA...done!\n")
		if(verbose) cat("RSDDA Lambda Parameter:", object$lambda, "\n")
	}
	else {
		object$lambda <- lambda
	}
	
	newdata <- data.matrix(newdata)
	
	if(!is.null(object$jointdiag.method) && object$jointdiag.method != "none") {
		newdata <- newdata %*% t(object$jointdiag.B)
	}
	
	predictions <- apply(newdata, 1, function(obs) {
		scores <- sapply(object$estimators, function(class_est) {
			var.rsdda <- (class_est$var.k)^(1 - object$lambda) * (class_est$var.pool)^(object$lambda)
			sum((obs - class_est$xbar)^2 * var.rsdda) - sum(log(var.rsdda)) - 2 * log(class_est$pi_k)
		})
		
		prediction <- object$classes[which.min(scores)]
		prediction
	})
	
	predictions <- factor(predictions, levels = object$classes)

	predictions
}