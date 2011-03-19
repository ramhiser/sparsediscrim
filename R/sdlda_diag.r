# Shrinkage-based Diagonal Linear Discriminant Analysis (SDLDA) with Covariance Matrix Diagonalization
# The DLDA classifier is a modification to LDA, where the off-diagonal elements
# of the pooled sample covariance matrix are set to zero.
# We use a shrinkage method based on Pang, Tong, and Zhao (2009).
# We assume the first column is named "labels" and holds a factor vector,
# which contains the class labels.
sdlda_diag <- function(train_df, num_alphas = 101, threshold = c("none", "hard"), tol = 0.01, ...) {
	obj <- list()
	
	threshold <- match.arg(threshold)
	
	N <- nrow(train_df)
	p <- ncol(train_df) - 1
	
	obj$N <- N
	obj$classes <- levels(train_df$labels)
	num_classes <- nlevels(train_df$labels)
	
	cov_pool <- Reduce('+', dlply(train_df, .(labels), function(df_k) {
		(n_k - 1) * cov(df_k[,-1])
	})) / N
	
	cov_eigen <- eigen(cov_pool, symmetric = TRUE)
	
	var_pool <- numeric(p)
	
	if(threshold == 'none') {
		obj$B <- cov_eigen$vectors
		var_pool <- cov_eigen$values
	} else if(threshold == 'hard') {
		q <- sum(cov_eigen$values >= tol)
		obj$B <- cov_eigen$vectors[, seq_len(q)]
		var_pool <- cov_eigen$values[seq_len(q)]
	}
	train_df <- cbind.data.frame(train_df$labels, data.matrix(train_df[, -1]) %*% obj$B)
	
	var_shrink <- var_shrinkage(N = N, K = num_classes, var.feature = var_pool, num_alphas = num_alphas, t = -1)
		
	estimators <- dlply(train_df, .(labels), function(df_k) {
		n_k <- nrow(df_k)
		pi_k <- n_k / N
		xbar <- as.vector(colMeans(df_k[,-1]))
		
		list(xbar = xbar, var = var_shrink, n_k = n_k, pi_k = pi_k)
	})

	class(obj) <- "sdlda_diag"
	
	obj
}

predict.sdlda_diag <- function(object, newdata) {
	if (!inherits(object, "sdlda_diag"))  {
		stop("object not of class 'sdlda_diag'")
	}
	newdata <- tcrossprod(data.matrix(newdata), object$B)
	
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