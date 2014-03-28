#' New Linear Discriminant Analysis CNLDA) Classifier
#'
#' TODO
#'
#' @export
#' @param x TODO
#' @param y TODO
#' @return obj TODO
nlda <- function(x, y) {
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
  Sig_eigen <- eigen(Sig, symmetric = TRUE)
  mean_eigenvals <- mean(Sig_eigen$values)
  Sig_eigen$values[Sig_eigen$values < mean_eigenvals] <- mean_eigenvals
  obj$Sig_inv <- Sig_eigen$vectors %*% diag(1/Sig_eigen$values) %*%
    t(Sig_eigen$vectors)

  class(obj) <- "nlda"

  return(obj)
}

#' Predicts the class membership of unlabeled observations with nlda classifier.
#'
#' # For a given nlda object, we predict the class of observations within the matrix newdata.
#' 
#' @export
#' @param object object of type 'nlda' that contains the trained classifier
#' @param newdata matrix containing the unlabeled observations to classify
#' @return list with predicted class and discriminant scores for each of the K classes
predict_nlda <- function(obj, newdata) {
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
