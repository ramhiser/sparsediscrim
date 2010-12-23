library(plyr)

# Shrinkage-based Diagonal Quadratic Discriminant Analysis (SDQDA)
# The SDQDA classifier is a modification to QDA, where the off-diagonal elements
# of the pooled sample covariance matrix are set to zero.
# We use a shrinkage method based on Pang, Tong, and Zhao (2009).
# We assume the first column is named "labels" and holds a factor vector,
# which contains the class labels.

# We assume the first column is named "labels" and holds a factor vector,
# which contains the class labels.
sdqda <- function(training.df, num.alphas = 5, jointdiag = "none", verbose = FALSE, ...) {
	sdqda.obj <- list()
	sdqda.obj$training <- training.df
	
	if(jointdiag != "none") {
		if(verbose) cat("Simultaneously diagonalizing covariance matrices\n")
		joint.diag.out <- joint.diagonalization(sdqda.obj$training, method = jointdiag)
		sdqda.obj$training <- joint.diag.out$transformed.df
		sdqda.obj$jointdiag.B <- joint.diag.out$B
		sdqda.obj$jointdiag.method <- joint.diag.out$method
		if(verbose) cat("Simultaneously diagonalizing covariance matrices...done!\n")
	}
	
	if(verbose) cat("Building SDQDA classifier\n")
	N <- nrow(sdqda.obj$training)
	num.classes <- nlevels(sdqda.obj$training$labels)
	
	estimators <- dlply(sdqda.obj$training, .(labels), function(class.df) {
		class.x <- as.matrix(class.df[, -1])
		dimnames(class.x) <- NULL
		
		n.k <- nrow(class.x)
		p.hat <- n.k / N
		
		xbar <- as.vector(colMeans(class.x))
		var <- apply(class.x, 2, function(col) {
			(n.k - 1) * var(col) / n.k
		})
		
		var.shrink <- var.shrinkage(N = N, K = 1, var.feature = var, num.alphas = num.alphas, t = -1)
		
		list(xbar = xbar, var = var.shrink, n = n.k, p.hat = p.hat)
	})
	if(verbose) cat("Building SDQDA classifier...done!\n")
	
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
	
	newdata <- as.matrix(newdata)
	dimnames(newdata) <- NULL
	
	predictions <- apply(newdata, 1, function(obs) {
		if(!is.null(object$jointdiag.method) && object$jointdiag.method != "none") {
			obs <- obs %*% t(object$jointdiag.B)
		}
		scores <- sapply(object$estimators, function(class.est) {
			sum((obs - class.est$xbar)^2 * class.est$var) - sum(log(class.est$var)) - 2 * log(class.est$p.hat)
		})
		predicted.class <- object$classes[which.min(scores)]
		predicted.class
	})
	
	predictions
}