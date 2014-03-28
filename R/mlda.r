#' Modified Linear Discriminant Analysis (MLDA) Classifier
#'
#' TODO
#'
#' @export
#' @reference TODO
#' @param x TODO
#' @param y TODO
#' @return obj TODO
mlda <- function(x, y) {
  obj <- list()
  y <- factor(y)
  obj$labels <- y
  obj$N <- length(y)
  obj$p <- ncol(x)
  obj$groups <- levels(y)
  obj$num_groups <- length(obj$groups)

  foreach_out <- foreach(k=levels(y)) %do% {
    est <- list()
    x_k <- x[which(y == k), ]
    est$label <- k
    est$n_k <- nrow(x_k)
    est$p_hat <- est$n_k / obj$N
    est$xbar_k <- colMeans(x_k)
    centered_x <- scale(x_k, center = TRUE, scale = FALSE)
    list(estimates = est, centered_x = centered_x)
  }
  obj$estimates <- lapply(foreach_out, function(k) k$estimates)
  centered_x <- do.call(rbind, lapply(foreach_out, function(k) k$centered_x))
  obj$Sig_inv <- invcov.shrink(centered_x, verbose = FALSE)
  class(obj) <- "mlda"
 
  return(obj)
}

#' Predicts the class membership of unlabeled observations with mlda classifier.
#'
#' # For a given mlda object, we predict the class of observations within the matrix newdata.
#' 
#' @export
#' @param object object of type 'mlda' that contains the trained classifier
#' @param newdata matrix containing the unlabeled observations to classify
#' @return list with predicted class and discriminant scores for each of the K classes
predict_mlda <- function(obj, newdata) {
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

