\name{rsdda}
\alias{rsdda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
TODO
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rsdda(train_df, num.alphas = 5, jointdiag = "none", verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{train_df}{
%%     ~~Describe \code{train_df} here~~
}
  \item{num.alphas}{
%%     ~~Describe \code{num.alphas} here~~
}
  \item{jointdiag}{
%%     ~~Describe \code{jointdiag} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(train_df, num.alphas = 5, jointdiag = "none", verbose = FALSE, ...) {
	rsdda.obj <- list()
	rsdda.obj$training <- train_df
	
	if(jointdiag != "none") {
		if(verbose) message("Simultaneously diagonalizing covariance matrices... ", appendLF = FALSE)
		joint.diag.out <- joint.diagonalization(rsdda.obj$training, method = jointdiag)
		rsdda.obj$training <- joint.diag.out$transformed.df
		rsdda.obj$jointdiag.B <- joint.diag.out$B
		rsdda.obj$jointdiag.method <- joint.diag.out$method
		if(verbose) message("done!")
	}

	if(verbose) cat("Building RSDDA classifier\n")
	N <- nrow(rsdda.obj$training)
	num.classes <- nlevels(rsdda.obj$training$labels)
	
	train_x <- as.matrix(rsdda.obj$training[,-1])
	dimnames(train_x) <- NULL
	
	estimators <- dlply(rsdda.obj$training, .(labels), function(df_k) {
		class.x <- as.matrix(df_k[, -1])
		dimnames(class.x) <- NULL
		
		n_k <- nrow(class.x)
		pi_k <- n_k / N
		
		xbar <- as.vector(colMeans(class.x))
		
		sum_squares <- apply(class.x, 2, function(col) {
			(n_k - 1) * var(col)
		})
		
		var <- apply(class.x, 2, function(col) {
			(n_k - 1) * var(col) / n_k
		})
		
		var.shrink <- var.shrinkage(N = n_k, K = 1, var.feature = var, num.alphas = num.alphas, t = -1)
		
		list(xbar = xbar, var.k = var.shrink, sum_squares = sum_squares, n = n_k, pi_k = pi_k)
	})
	
	var_pool <- colSums(laply(estimators, function(class_est) class_est$sum_squares)) / N
	var.pool.shrink <- var.shrinkage(N = N, K = num.classes, var.feature = var_pool, num.alphas = num.alphas, t = -1)
	
	estimators <- llply(estimators, function(class_estimators) {
		class_estimators$var.pool <- var.pool.shrink
		class_estimators
	})
	
	if(verbose) cat("Building RSDDA classifier...done!\n")
	
	rsdda.obj$N <- N
	rsdda.obj$classes <- levels(rsdda.obj$training$labels)
	rsdda.obj$estimators <- estimators
	
	class(rsdda.obj) <- "rsdda"

	rsdda.obj
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
