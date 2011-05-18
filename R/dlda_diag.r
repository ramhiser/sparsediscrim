# Diagonalized Linear Discriminant Analysis (DLDA) with Covariance Matrix Diagonalization
# The DLDA classifier is a modification to LDA, where the off-diagonal elements
# of the pooled sample covariance matrix are set to zero.
# We assume the first column is named "labels" and holds a factor vector,
# which contains the class labels.
dlda_diag <- function(train_df, threshold = c("none", "eigen_pct", "hard", "soft"), delta = 0.01, eigen_pct = 0.95, ...) {
	obj <- list()
	
	threshold <- match.arg(threshold)
	obj$threshold <- threshold
	
	N <- nrow(train_df)
	p <- ncol(train_df) - 1
	obj$N <- N
	obj$classes <- levels(train_df$labels)
	
	cov_pool <- Reduce('+', dlply(train_df, .(labels), function(df_k) {
		(nrow(df_k) - 1) * cov(df_k[,-1])
	})) / N
	
	cov_eigen <- eigen(cov_pool, symmetric = TRUE)
	var_pool <- numeric(p)
	
	if(threshold == 'none') {
		obj$B <- cov_eigen$vectors
		var_pool <- cov_eigen$values
	} else if(threshold == 'hard') {
		q <- sum(cov_eigen$values >= delta)
		obj$B <- cov_eigen$vectors[, seq_len(q)]
		var_pool <- cov_eigen$values[seq_len(q)]
	} else if(threshold == 'eigen_pct') {
		q <- sum(cumsum(cov_eigen$values) / sum(cov_eigen$value) <= eigen_pct)
		obj$B <- cov_eigen$vectors[, seq_len(q)]
		var_pool <- cov_eigen$values[seq_len(q)]
	} else if(threshold == 'soft') {
		shrunken_eigenvals <- cov_eigen$values - delta
		shrunken_eigenvals[which(shrunken_eigenvals < 0)] <- 0
		q <- sum(shrunken_eigenvals > 0)
		obj$B <- cov_eigen$vectors[, seq_len(q)]
		var_pool <- shrunken_eigenvals[seq_len(q)]
	}
	train_df <- cbind.data.frame(labels = train_df$labels, data.matrix(train_df[, -1]) %*% obj$B)

	obj$estimators <- dlply(train_df, .(labels), function(df_k) {
		n_k <- nrow(df_k)
		pi_k <- n_k / N
		if(is.vector(df_k[,-1])) {
			xbar <- mean(df_k[,-1])
		} else {
			xbar <- as.vector(colMeans(df_k[,-1]))
		}
		list(xbar = xbar, var = var_pool, n_k = n_k, pi_k = pi_k)
	})

	class(obj) <- "dlda_diag"
	
	obj
}

dlda_diag.cv <- function(train_df, threshold = c("soft", "hard"), hold_out = 5, delta = seq(0.1, 0.5, by = 0.1), verbose = FALSE) {
	laply(cv_partition(train_df$labels, k = hold_out), function(held_out) {
		cv_train_df <- train_df[-held_out,]
		cv_test_df <- train_df[held_out,]
		
		errors <- sapply(delta, function(del) {
			cv_out <- dlda_diag(cv_train_df, threshold = threshold, delta = del)
			cv_predictions <- predict.dlda_diag(cv_out, cv_test_df[,-1])
			sum(cv_test_df$labels != cv_predictions)
		})
		errors
	})
}

cv_out <- dlda_diag(cv_train_df, threshold = threshold, delta = 0.2)
cv_predictions <- predict.dlda_diag(cv_out, cv_test_df[,-1])
sum(cv_test_df$labels != cv_predictions)

predict.dlda_diag <- function(object, newdata) {
	if (!inherits(object, "dlda_diag"))  {
		stop("object not of class 'dlda_diag'")
	}
	if(is.vector(newdata)) {
		newdata <- matrix(data.matrix(newdata), nrow = 2)
	} else {
		newdata <- data.matrix(newdata)
	}
	newdata <- newdata %*% object$B

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