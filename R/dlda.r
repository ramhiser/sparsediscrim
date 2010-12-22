library(plyr)

# Diagonalized Linear Discriminant Analysis (DLDA)
# The DLDA classifier is a modification to LDA, where the off-diagonal elements
# of the pooled sample covariance matrix are set to zero.
# We assume the first column is named "labels" and holds a factor vector,
# which contains the class labels.
dlda <- function(training.df, jointdiag = "none") {
	N <- nrow(training.df)
	
	training.x <- as.matrix(training.df[,-1])
	dimnames(training.x) <- NULL
	pooled.var <- apply(training.x, 2, function(col) {
		(N - 1) * var(col) / N
	})
	
	estimators <- dlply(training.df, .(labels), function(class.df) {
		class.x <- as.matrix(class.df[,-1])
		dimnames(class.x) <- NULL
		
		n.k <- nrow(class.df)
		p.hat <- n.k / N
		xbar <- as.vector(colMeans(class.x))
		var <- pooled.var
		list(xbar = xbar, var = var, n = n.k, p.hat = p.hat)
	})
	
	list(N = N, classes = levels(training.df$labels), estimators = estimators)
}

predict.dlda <- function(object, newdata) {
	newdata <- as.matrix(newdata)
	dimnames(newdata) <- NULL
	
	predictions <- apply(newdata, 1, function(obs) {
		scores <- sapply(object$estimators, function(class.est) {
			sum((obs - class.est$xbar)^2 / class.est$var) - 2 * log(class.est$p.hat)
		})

		predicted.class <- object$classes[which.min(scores)]
		predicted.class
	})
	
	predictions
}