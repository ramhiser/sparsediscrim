#' Linear Discriminant Analysis with the Moore-Penrose Pseudo-Inverse
#'
#' When the pooled sample covariance matrix is singular, the linear discriminant
#' function is incalculable. A common method to overcome this issue is to
#' replace the inverse of the pooled sample covariance matrix with the
#' Moore-Penrose (MP) pseudo-inverse, which is unique and always exists, and
#' furthermore, when the pooled sample covariance matrix is nonsingular, it
#' is equal to the MP pseudo-inverse.
#'
#' @export
#' @param x TODO
#' @param y TODO
#' @return obj TODO
lda_pseudo <- function(x, y) {
  obj <- list()
  y <- factor(y)
  obj$labels <- y
  obj$N <- length(y)
  obj$p <- ncol(x)
  obj$groups <- levels(y)
  obj$num_groups <- length(obj$groups)

  obj$estimates <- foreach(k=levels(y)) %do% {
    est <- list()
    x_k <- x[which(y == k), ]
    est$label <- k
    est$n_k <- nrow(x_k)
    est$p_hat <- est$n_k / obj$N
    est$xbar_k <- colMeans(x_k)
    est
  }
  Sig <- foreach(k=levels(y)) %do% {
    x_k <- x[which(y == k), ]
    (nrow(x_k) - 1) * cov(x_k)
  }
  Sig <- Reduce('+', Sig) / (obj$N - obj$num_groups)
  obj$Sig_inv <- try(solve(Sig), silent = TRUE)
  if(!is.matrix(obj$Sig_inv)) obj$Sig_inv <- ginv(Sig)

  class(obj) <- "lda_pseudo"

  return(obj)
}

#' Predicts the class membership of unlabeled observations with lda_pseudo classifier.
#'
#' # For a given lda_pseudo object, we predict the class of observations within the matrix newdata.
#' 
#' @export
#' @param object object of type 'lda_pseudo' that contains the trained classifier
#' @param newdata matrix containing the unlabeled observations to classify
#' @return list with predicted class and discriminant scores for each of the K classes
predict_lda_pseudo <- function(obj, newdata) {
  if(is.vector(newdata)) newdata <- matrix(newdata, nrow = 1)
  scores <- sapply(obj$estimates, function(est) {
    apply(newdata, 1, function(obs) {
      anderson_ldf(obs, est$xbar_k, obj$Sig_inv, est$p_hat)
    })
  })
	
  if(is.vector(scores)) {
    min_scores <- which.min(scores)
  } else {
    min_scores <- apply(scores, 1, which.min)
  }

  predicted <- factor(obj$groups[min_scores], levels = obj$groups)

  list(scores = scores, predicted = predicted)
}
