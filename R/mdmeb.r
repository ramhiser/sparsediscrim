# Minimum Distance Rule using Modified Empirical Bayes (MDMEB) from Srivistava and Kubokawa (2007)
#'
#' TODO
#'
#' @export
#' @reference TODO
#' @param x TODO
#' @param y TODO
#' @return obj TODO
mdmeb <- function(x, y, eigen_pct = 0.95) {
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

  pctile <- quantile(Sig_eigen$values, probs = 1 - eigen_pct)
  kept_evals <- which(Sig_eigen$values >= pctile) 
  # TODO: Are these e'vals shrunken correctly?
  shrunk_evals <- Sig_eigen$values[kept_evals] + mean(Sig_eigen$values[kept_evals])

  obj$Sig_inv <- tcrossprod(Sig_eigen$vectors[,kept_evals] %*% diag(1/shrunk_evals), Sig_eigen$vectors[,kept_evals])
  class(obj) <- "mdmeb"

  return(obj)
}

#' Predicts the class membership of unlabeled observations with mdmeb classifier.
#'
#' # For a given mdmeb object, we predict the class of observations within the matrix newdata.
#' 
#' @export
#' @param object object of type 'mdmeb' that contains the trained classifier
#' @param newdata matrix containing the unlabeled observations to classify
#' @return list with predicted class and discriminant scores for each of the K classes
predict_mdmeb <- function(obj, newdata) {
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

