# Diagonalized Linear Discriminant Analysis (DLDA)
# The DLDA classifier is a modification to LDA, where the off-diagonal elements
# of the pooled sample covariance matrix are set to zero.
# We assume the first column is named "labels" and holds a factor vector,
# which contains the class labels.
dlda <- function(training.df, jointdiag = "none", verbose = FALSE, ...) {
	dlda.obj <- list()
	dlda.obj$training <- training.df
	
	if(jointdiag != "none") {
		if(verbose) cat("Simultaneously diagonalizing covariance matrices\n")
		joint.diag.out <- joint.diagonalization(dlda.obj$training, method = jointdiag)
		dlda.obj$training <- joint.diag.out$transformed.df
		dlda.obj$jointdiag.B <- joint.diag.out$B
		dlda.obj$jointdiag.method <- joint.diag.out$method
		if(verbose) cat("Simultaneously diagonalizing covariance matrices...done!\n")
	}

	if(verbose) cat("Building DLDA classifier\n")
	N <- nrow(dlda.obj$training)
	
	training.x <- as.matrix(dlda.obj$training[,-1])
	dimnames(training.x) <- NULL
	
	estimators <- dlply(dlda.obj$training, .(labels), function(class.df) {
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
	
	estimators <- llply(estimators, function(class.estimators) {
		class.estimators$var <- var.pooled
		class.estimators
	})
	if(verbose) cat("Building DLDA classifier...done!\n")
	
	dlda.obj$N <- N
	dlda.obj$classes <- levels(dlda.obj$training$labels)
	dlda.obj$estimators <- estimators
	
	class(dlda.obj) <- "dlda"
	
	dlda.obj
}

predict.dlda <- function(object, newdata) {
	if (!inherits(object, "dlda"))  {
		stop("object not of class 'dlda'")
	}
	newdata <- as.matrix(newdata)
	dimnames(newdata) <- NULL
	
	predictions <- apply(newdata, 1, function(obs) {
		if(!is.null(object$jointdiag.method) && object$jointdiag.method != "none") {
			obs <- obs %*% t(object$jointdiag.B)
		}
		scores <- sapply(object$estimators, function(class.est) {
			sum((obs - class.est$xbar)^2 / class.est$var) - 2 * log(class.est$p.hat)
		})
		predicted.class <- object$classes[which.min(scores)]
		predicted.class
	})
	
	predictions
}