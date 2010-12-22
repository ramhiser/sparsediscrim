library(plyr)

# Regularized Shrinkage-based Diagonal Discriminant Analysis (RSDDA)
# The RSDDA classifier is a regularized version of SDQDA, where the off-diagonal elements
# of the pooled sample covariance matrix are set to zero.
# We use a shrinkage method based on Pang, Tong, and Zhao (2009).
# We assume the first column is named "labels" and holds a factor vector,
# which contains the class labels.

# We assume the first column is named "labels" and holds a factor vector,
# which contains the class labels.

# num.alphas: the grid size of alphas considered in estimating each shrinkage parameter alpha
rsdda <- function(training.df, num.alphas = 5, jointdiag = "none") {
	N <- nrow(training.df)
	num.classes <- nlevels(training.df$labels)
	
	training.x <- as.matrix(training.df[,-1])
	dimnames(training.x) <- NULL
	var.pooled <- apply(training.x, 2, function(col) {
		(N - 1) * var(col) / N
	})
	
	var.pool.shrink <- var.shrinkage(N = N, K = num.classes, var.feature = var.pooled, num.alphas = num.alphas, t = -1)
	
	estimators <- dlply(training.df, .(labels), function(class.df) {
		x.k <- as.matrix(class.df[, -1])
		dimnames(x.k) <- NULL
		
		n.k <- nrow(x.k)
		p.hat <- n.k / N
		
		xbar <- as.vector(colMeans(x.k))
		var <- apply(x.k, 2, function(col) {
			(n.k - 1) * var(col) / n.k
		})
		
		var.shrink <- var.shrinkage(N = N, K = 1, var.feature = var, num.alphas = num.alphas, t = -1)
		
		list(xbar = xbar, var.k = var.shrink, var.pool = var.pool.shrink, n = n.k, p.hat = p.hat)
	})
	
	list(N = N, classes = levels(training.df$labels), estimators = estimators)
}

# Estimates the optimal value of lambda used in RSDDA by leave-k-out crossvalidation.
# grid.size: the grid size of lambdas considered in estimating the regularization parameter lambda
# k: Leave-k-out crossvalidation is used to estimate the optimal value of lambda.
model.select.rsdda <- function(training.df, object, grid.size = 5, k = 1) {
	# Estimate lambda with Leave-k-out crossvalidation.
	if(k > 1) {
		k <- 1
		warning("Leave-k-out crossvalidation for k > 1 is not implemented yet. Using k = 1 instead...\n")
	}
	
	lambdas <- seq(0, 1, length = grid.size)
	
	# Calculate the Leave-K-out Error Rates for each lambda in the grid.
	predictions <- laply(seq_len(object$N), function(i) {
		loo.rsdda.obj <- rsdda(training.df[-i, ], num.alphas = grid.size)
		laply(lambdas, function(lambda) {
			loo.rsdda.obj$lambda <- lambda
			predict.rsdda(loo.rsdda.obj, training.df[i, -1])
		})
	})

	error.rates.loo <- apply(predictions, 2, function(predictions.lambda) {
		mean(training.df[, 1] != predictions.lambda)
	})
	
	# Find the lambdas that attain the minimum Leave-K-out error rate.
	optimal.lambdas <- lambdas[which(error.rates.loo == min(error.rates.loo))]
	
	# Finally, we add the optimal lambda to the RSDDA object.
	# If there is a tie for optimal lambda, randomly choose one from the optimal ones.
	object$lambda <- sample(optimal.lambdas, 1)
	
	return(object)
}

predict.rsdda <- function(object, newdata) {
	newdata <- as.matrix(newdata)
	dimnames(newdata) <- NULL
	
	predictions <- apply(newdata, 1, function(obs) {
		scores <- sapply(object$estimators, function(class.est) {
			var.rsdda <- (class.est$var.k)^(1 - object$lambda) * (class.est$var.pool)^(object$lambda)
			sum((obs - class.est$xbar)^2 * var.rsdda) - sum(log(var.rsdda)) - 2 * log(class.est$p.hat)
		})
		
		predicted.class <- object$classes[which.min(scores)]
		predicted.class
	})

	predictions
}