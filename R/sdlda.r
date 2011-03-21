# Shrinkage-based Diagonal Linear Discriminant Analysis (SDLDA)
# The DLDA classifier is a modification to LDA, where the off-diagonal elements
# of the pooled sample covariance matrix are set to zero.
# We use a shrinkage method based on Pang, Tong, and Zhao (2009).
# We assume the first column is named "labels" and holds a factor vector,
# which contains the class labels.

# We assume the first column is named "labels" and holds a factor vector,
# which contains the class labels.
sdlda <- function(train_df, num_alphas = 101) {
	obj <- list()
	
	N <- nrow(train_df)
	obj$training <- train_df
	obj$N <- N
	obj$classes <- levels(train_df$labels)
	obj$num_classes <- nlevels(train_df$labels)
	
	estimators <- dlply(obj$training, .(labels), function(df_k) {
		n_k <- nrow(df_k)
		pi_k <- n_k / N
		xbar <- as.vector(colMeans(df_k[,-1]))
		sum_squares <- (n_k - 1) * apply(df_k[,-1], 2, var)
		list(xbar = xbar, sum_squares = sum_squares, n_k = n_k, pi_k = pi_k)
	})
	
	# TODO: Calculate var_pool before calculating other estimators,
	#		so that sum_squares is not being carried around for each class.
	var_pool <- colSums(laply(estimators, function(class_est) class_est$sum_squares)) / N
	var_shrink <- var_shrinkage(N = N, K = num_classes, var_feature = var_pool, num_alphas = num_alphas, t = -1)
	
	obj$estimators <- llply(estimators, function(class_estimators) {
		class_estimators$var <- var_shrink
		class_estimators
	})

	class(obj) <- "sdlda"
	
	obj
}

predict.sdlda <- function(object, newdata) {
	if (!inherits(object, "sdlda"))  {
		stop("object not of class 'sdlda'")
	}
	
	newdata <- data.matrix(newdata)

	predictions <- apply(newdata, 1, function(obs) {
		scores <- sapply(object$estimators, function(class_est) {
			sum((obs - class_est$xbar)^2 * class_est$var) - 2 * log(class_est$pi_k)
		})
		prediction <- object$classes[which.min(scores)]
		prediction
	})
	
	predictions <- factor(predictions, levels = object$classes)
	
	predictions
}