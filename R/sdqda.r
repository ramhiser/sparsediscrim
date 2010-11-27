library(plyr)

# Shrinkage-based Diagonal Quadratic Discriminant Analysis (SDQDA)
# The SDQDA classifier is a modification to QDA, where the off-diagonal elements
# of the pooled sample covariance matrix are set to zero.
# We use a shrinkage method based on Pang, Tong, and Zhao (2009).
# We assume the first column is named "labels" and holds a factor vector,
# which contains the class labels.

# We assume the first column is named "labels" and holds a factor vector,
# which contains the class labels.
sdqda <- function(training.df, num.alphas = 5) {
	N <- nrow(training.df)
	num.classes <- nlevels(training.df$labels)
	
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
		
		list(xbar = xbar, var = var.shrink, n = n.k, p.hat = p.hat)
	})
	
	list(N = N, classes = levels(training.df$labels), estimators = estimators)
}

predict.sdqda <- function(object, newdata) {
	newdata <- as.matrix(newdata)
	dimnames(newdata) <- NULL
	
	predictions <- apply(newdata, 1, function(obs) {
		scores <- laply(object$estimators, function(class.est) {
			sum((obs - class.est$xbar)^2 * class.est$var) - sum(log(class.est$var)) - 2 * log(class.est$p.hat)
		})
		predicted.class <- object$classes[which.min(scores)]
		predicted.class
	})
	
	predictions
}