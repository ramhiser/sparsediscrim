#' Calculates the percentage change for each element in the given vector
#'
#' For each element in the vector, \code{x}, we calculate the percentage
#' change.
#'
#' The returned vector has one less observation than \code{x}.
#' @param x vector
#' @return vector of length \code{length(x) - 1}
#' @examples
#' x <- seq_len(10)^2
#' out <- pct_change(x)
#' out
#' length(out) == 9
pct_change <- function(x) {
  x <- as.vector(x)
  diff(x) / head(x, length(x) - 1)
}
