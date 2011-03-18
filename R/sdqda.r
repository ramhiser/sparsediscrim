# Shrinkage-based Diagonal Quadratic Discriminant Analysis (SDQDA)
# The SDQDA classifier is a modification to QDA, where the off-diagonal elements
# of the pooled sample covariance matrix are set to zero.
# We use a shrinkage method based on Pang, Tong, and Zhao (2009).
# We assume the first column is named "labels" and holds a factor vector,
# which contains the class labels.

# We assume the first column is named "labels" and holds a factor vector,
# which contains the class labels.
sdqda <- function(train_df, num.alphas = 5, jointdiag = "none", verbose = FALSE, ...) {
	sdqda.obj <- list()
	sdqda.obj$training <- train_df
	
	if(jointdiag != "none") {
		if(verbose) message("Simultaneously diagonalizing covariance matrices... ", appendLF = FALSE)
		joint.diag.out <- joint.diagonalization(sdqda.obj$training, method = jointdiag)
		sdqda.obj$training <- joint.diag.out$transformed.df
		sdqda.obj$jointdiag.B <- joint.diag.out$B
		sdqda.obj$jointdiag.method <- joint.diag.out$method
		if(verbose) message("done!")
	}
	
	if(verbose) message("Building SDQDA classifier... ", appendLF = FALSE)
	N <- nrow(sdqda.obj$training)
	num.classes <- nlevels(sdqda.obj$training$labels)
	
	estimators <- dlply(sdqda.obj$training, .(labels), function(df_k) {
		class.x <- as.matrix(df_k[, -1])
		dimnames(class.x) <- NULL
		
		n_k <- nrow(class.x)
		pi_k <- n_k / N
		
		xbar <- as.vector(colMeans(class.x))
		var <- apply(class.x, 2, function(col) {
			(n_k - 1) * var(col) / n_k
		})
		
		var.shrink <- var.shrinkage(N = n_k, K = 1, var.feature = var, num.alphas = num.alphas, t = -1)
		
		list(xbar = xbar, var = var.shrink, n = n_k, pi_k = pi_k)
	})
	if(verbose) message("done!")
	
	sdqda.obj$N <- N
	sdqda.obj$classes <- levels(sdqda.obj$training$labels)
	sdqda.obj$estimators <- estimators
	
	class(sdqda.obj) <- "sdqda"
	
	sdqda.obj
	
}

predict.sdqda <- function(object, newdata) {
	if (!inherits(object, "sdqda"))  {
		stop("object not of class 'sdqda'")
	}
	
	newdata <- data.matrix(newdata)
	
	if(!is.null(object$jointdiag.method) && object$jointdiag.method != "none") {
		newdata <- newdata %*% t(object$jointdiag.B)
	}
	
	predictions <- apply(newdata, 1, function(obs) {
		scores <- sapply(object$estimators, function(class_est) {
			sum((obs - class_est$xbar)^2 * class_est$var) - sum(log(class_est$var)) - 2 * log(class_est$pi_k)
		})
		predicted.class <- object$classes[which.min(scores)]
		predicted.class
	})
	
	predictions
}