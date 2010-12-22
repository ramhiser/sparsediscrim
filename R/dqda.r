library(plyr)

# Diagonalized Quadratic Discriminant Analysis (DQDA)
# The DQDA classifier is a modification to LDA, where the off-diagonal elements
# of the pooled sample covariance matrix are set to zero.
# We assume the first column is named "labels" and holds a factor vector,
# which contains the class labels.

# We assume the first column is named "labels" and holds a factor vector,
# which contains the class labels.
dqda <- function(training.df, jointdiag = "none") {
	N <- nrow(training.df)
	
	estimators <- dlply(training.df, .(labels), function(class.df) {
		n.k <- nrow(class.df)
		p.hat <- n.k / N
		xbar <- as.vector(colMeans(class.df[, -1]))
		var <- apply(class.df[,-1], 2, function(col) {
			(n.k - 1) * var(col) / n.k
		})
		list(xbar = xbar, var = var, n = n.k, p.hat = p.hat)
	})
	
	list(N = N, classes = levels(training.df$labels), estimators = estimators)
}

predict.dqda <- function(object, newdata) {
	newdata <- as.matrix(newdata)
	dimnames(newdata) <- NULL
	
	predictions <- apply(newdata, 1, function(obs) {
		scores <- sapply(object$estimators, function(class.est) {
			sum((obs - class.est$xbar)^2 / class.est$var) + sum(log(class.est$var)) - 2 * log(class.est$p.hat)
		})
		predicted.class <- object$classes[which.min(scores)]
		predicted.class
	})
	
	predictions
}