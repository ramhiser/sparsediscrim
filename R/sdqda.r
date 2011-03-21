# Shrinkage-based Diagonal Quadratic Discriminant Analysis (SDQDA)
# The SDQDA classifier is a modification to QDA, where the off-diagonal elements
# of the pooled sample covariance matrix are set to zero.
# We use a shrinkage method based on Pang, Tong, and Zhao (2009).
# We assume the first column is named "labels" and holds a factor vector,
# which contains the class labels.

# We assume the first column is named "labels" and holds a factor vector,
# which contains the class labels.
sdqda <- function(train_df, num_alphas = 101) {
	obj <- list()
	
	N <- nrow(train_df)
	obj$training <- train_df
	obj$N <- N
	obj$classes <- levels(train_df$labels)
	obj$num_classes <- nlevels(train_df$labels)
	
	obj$estimators <- dlply(obj$training, .(labels), function(df_k) {
		n_k <- nrow(df_k)
		pi_k <- n_k / N
		xbar <- as.vector(colMeans(df_k[,-1]))
		var <- (n_k - 1) / n_k * apply(df_k[,-1], 2, var)
		var_shrink <- var_shrinkage(N = n_k, K = 1, var_feature = var, num_alphas = num_alphas, t = -1)
		list(xbar = xbar, var = var_shrink, n_k = n_k, pi_k = pi_k)
	})

	class(obj) <- "sdqda"
	
	obj	
}

predict.sdqda <- function(object, newdata) {
	if (!inherits(object, "sdqda"))  {
		stop("object not of class 'sdqda'")
	}
	
	if(is.vector(newdata)) {
		newdata <- matrix(data.matrix(newdata), nrow = 2)
	} else {
		newdata <- data.matrix(newdata)
	}
	
	predictions <- apply(newdata, 1, function(obs) {
		scores <- sapply(object$estimators, function(class_est) {
			sum((obs - class_est$xbar)^2 * class_est$var) - sum(log(class_est$var)) - 2 * log(class_est$pi_k)
		})
		prediction <- object$classes[which.min(scores)]
		prediction
	})
	
	predictions <- factor(predictions, levels = object$classes)
	
	predictions
}