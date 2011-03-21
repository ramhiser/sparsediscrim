# Diagonalized Quadratic Discriminant Analysis (DQDA)
# The DQDA classifier is a modification to LDA, where the off-diagonal elements
# of the pooled sample covariance matrix are set to zero.
# We assume the first column is named "labels" and holds a factor vector,
# which contains the class labels.

# We assume the first column is named "labels" and holds a factor vector,
# which contains the class labels.
dqda <- function(train_df, jointdiag = "none", verbose = FALSE, ...) {
	obj <- list()
	
	N <- nrow(train_df)
	obj$training <- train_df
	obj$N <- N
	obj$classes <- levels(train_df$labels)
	
	obj$estimators <- dlply(obj$training, .(labels), function(df_k) {
		n_k <- nrow(df_k)
		pi_k <- n_k / N
		xbar <- as.vector(colMeans(df_k[,-1]))
		var <- (n_k - 1) / n_k * apply(df_k[,-1], 2, var)
		list(xbar = xbar, var = var, n_k = n_k, pi_k = pi_k)
	})

	class(obj) <- "dqda"
	
	obj
}

predict.dqda <- function(object, newdata) {
	if (!inherits(object, "dqda"))  {
		stop("object not of class 'dqda'")
	}
	newdata <- data.matrix(newdata)

	predictions <- apply(newdata, 1, function(obs) {
		scores <- sapply(object$estimators, function(class_est) {
			sum((obs - class_est$xbar)^2 / class_est$var) + sum(log(class_est$var)) - 2 * log(class_est$pi_k)
		})
		prediction <- object$classes[which.min(scores)]
		prediction
	})
	
	predictions <- factor(predictions, levels = object$classes)
	
	predictions
}