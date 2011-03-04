# Shrinkage-based Diagonal Linear Discriminant Analysis (SDLDA)
# The DLDA classifier is a modification to LDA, where the off-diagonal elements
# of the pooled sample covariance matrix are set to zero.
# We use a shrinkage method based on Pang, Tong, and Zhao (2009).
# We assume the first column is named "labels" and holds a factor vector,
# which contains the class labels.

# We assume the first column is named "labels" and holds a factor vector,
# which contains the class labels.
sdlda <- function(training.df, num.alphas = 5, jointdiag = "none", verbose = FALSE, ...) {
	sdlda.obj <- list()
	sdlda.obj$training <- training.df
	
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
	
	training.x <- as.matrix(sdlda.obj$training[,-1])
	dimnames(training.x) <- NULL
	
	estimators <- dlply(sdlda.obj$training, .(labels), function(class.df) {
		class.x <- as.matrix(class.df[,-1])
		dimnames(class.x) <- NULL
		
		n.k <- nrow(class.df)
		p.hat <- n.k / N
		xbar <- as.vector(colMeans(class.x))
		
		sum.squares <- apply(class.x, 2, function(col) {
			(n.k - 1) * var(col)
		})
		
		list(xbar = xbar, sum.squares = sum.squares, n = n.k, p.hat = p.hat)
	})
	
	var.pooled <- colSums(laply(estimators, function(class.est) class.est$sum.squares)) / N
	var.shrink <- var.shrinkage(N = N, K = num.classes, var.feature = var.pooled, num.alphas = num.alphas, t = -1)
	
	estimators <- llply(estimators, function(class.estimators) {
		class.estimators$var <- var.shrink
		class.estimators
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
	
	newdata <- as.matrix(newdata)
	dimnames(newdata) <- NULL
	
	predictions <- apply(newdata, 1, function(obs) {
		if(!is.null(object$jointdiag.method) && object$jointdiag.method != "none") {
			obs <- obs %*% t(object$jointdiag.B)
		}
		scores <- sapply(object$estimators, function(class.est) {
			sum((obs - class.est$xbar)^2 * class.est$var) - 2 * log(class.est$p.hat)
		})
		predicted.class <- object$classes[which.min(scores)]
		predicted.class
	})
	
	predictions
}