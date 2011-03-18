# Shrinkage-based Diagonal Linear Discriminant Analysis (SDLDA)
# The DLDA classifier is a modification to LDA, where the off-diagonal elements
# of the pooled sample covariance matrix are set to zero.
# We use a shrinkage method based on Pang, Tong, and Zhao (2009).
# We assume the first column is named "labels" and holds a factor vector,
# which contains the class labels.

# We assume the first column is named "labels" and holds a factor vector,
# which contains the class labels.
sdlda <- function(train_df, num.alphas = 5, jointdiag = "none", verbose = FALSE, ...) {
	sdlda.obj <- list()
	sdlda.obj$training <- train_df
	
	if(jointdiag != "none") {
		if(verbose) message("Simultaneously diagonalizing covariance matrices... ", appendLF = FALSE)
		joint.diag.out <- joint.diagonalization(sdlda.obj$training, method = jointdiag)
		sdlda.obj$training <- joint.diag.out$transformed.df
		sdlda.obj$jointdiag.B <- joint.diag.out$B
		sdlda.obj$jointdiag.method <- joint.diag.out$method
		if(verbose) message("done!")
	}
	
	if(verbose) message("Building SDLDA classifier... ", appendLF = FALSE)
	N <- nrow(sdlda.obj$training)
	num.classes <- nlevels(sdlda.obj$training$labels)
	
	train_x <- as.matrix(sdlda.obj$training[,-1])
	dimnames(train_x) <- NULL
	
	estimators <- dlply(sdlda.obj$training, .(labels), function(df_k) {
		class.x <- as.matrix(df_k[,-1])
		dimnames(class.x) <- NULL
		
		n_k <- nrow(df_k)
		pi_k <- n_k / N
		xbar <- as.vector(colMeans(class.x))
		
		sum_squares <- apply(class.x, 2, function(col) {
			(n_k - 1) * var(col)
		})
		
		list(xbar = xbar, sum_squares = sum_squares, n = n_k, pi_k = pi_k)
	})
	
	var_pool <- colSums(laply(estimators, function(class_est) class_est$sum_squares)) / N
	var.shrink <- var.shrinkage(N = N, K = num.classes, var.feature = var_pool, num.alphas = num.alphas, t = -1)
	
	estimators <- llply(estimators, function(class_estimators) {
		class_estimators$var <- var.shrink
		class_estimators
	})
	if(verbose) message("done!")
	
	sdlda.obj$N <- N
	sdlda.obj$classes <- levels(sdlda.obj$training$labels)
	sdlda.obj$estimators <- estimators
	
	class(sdlda.obj) <- "sdlda"
	
	sdlda.obj
}

predict.sdlda <- function(object, newdata) {
	if (!inherits(object, "sdlda"))  {
		stop("object not of class 'sdlda'")
	}
	
	newdata <- data.matrix(newdata)
	
	if(!is.null(object$jointdiag.method) && object$jointdiag.method != "none") {
		newdata <- newdata %*% t(object$jointdiag.B)
	}
	
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