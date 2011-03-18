# Diagonalized Linear Discriminant Analysis (DLDA)
# The DLDA classifier is a modification to LDA, where the off-diagonal elements
# of the pooled sample covariance matrix are set to zero.
# We assume the first column is named "labels" and holds a factor vector,
# which contains the class labels.
dlda <- function(train_df) {
	obj <- list()
	
	N <- nrow(train_df)
	
	obj$training <- train_df
	obj$N <- N
	obj$classes <- levels(train_df$labels)
	
	estimators <- dlply(obj$training, .(labels), function(df_k) {
		n_k <- nrow(df_k)
		pi_k <- n_k / N
		xbar <- as.vector(colMeans(df_k[,-1]))
		
		sum_squares <- apply(x_k, 2, function(col) {
			(n_k - 1) * var(col)
		})
		
		list(xbar = xbar, sum_squares = sum_squares, n_k = n_k, pi_k = pi_k)
	})
	
	var_pool <- colSums(laply(estimators, function(class_est) class_est$sum_squares)) / N
	
	obj$estimators <- llply(estimators, function(class_estimators) {
		class_estimators$var <- var_pool
		class_estimators
	})

	class(obj) <- "dlda"
	
	obj
}

predict.dlda <- function(object, newdata) {
	if (!inherits(object, "dlda"))  {
		stop("object not of class 'dlda'")
	}
	newdata <- data.matrix(newdata)
	
	predictions <- apply(newdata, 1, function(obs) {
		scores <- sapply(object$estimators, function(class_est) {
			sum((obs - class_est$xbar)^2 / class_est$var) - 2 * log(class_est$pi_k)
		})
		prediction <- object$classes[which.min(scores)]
		prediction
	})
	
	predictions <- factor(predictions, levels = object$classes)
	
	predictions
}