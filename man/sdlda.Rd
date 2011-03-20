\name{sdlda}
\alias{sdlda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
TODO
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sdlda(train_df, num_alphas = 5, jointdiag = "none", verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{train_df}{
%%     ~~Describe \code{train_df} here~~
}
  \item{num_alphas}{
%%     ~~Describe \code{num_alphas} here~~
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
function(train_df, num_alphas = 5, jointdiag = "none", verbose = FALSE, ...) {
	sdlda.obj <- list()
	sdlda.obj$training <- train_df
	
	if(jointdiag != "none") {
		if(verbose) message("Simultaneously diagonalizing covariance matrices... ", appendLF = FALSE)
		joint.diag.out <- joint.diagonalization(sdlda.obj$training, method = jointdiag)
		sdlda.obj$training <- joint.diag.out$transformed.df
		sdlda.obj$jointdiag.B <- joint.diag.out$B
		sdlda.obj$jointdiag.method <- joint.diag.out$method
		if(verbose) message("done!")
	}
	
	if(verbose) cat("Building SDLDA classifier\n")
	N <- nrow(sdlda.obj$training)
	num.classes <- nlevels(sdlda.obj$training$labels)
	
	train_x <- as.matrix(sdlda.obj$training[,-1])
	dimnames(train_x) <- NULL
	
	estimators <- dlply(sdlda.obj$training, .(labels), function(df_k) {
		class.x <- as.matrix(df_k[,-1])
		dimnames(class.x) <- NULL
		
		n_k <- nrow(df_k)
		pi_k <- n_k / N
		xbar <- as.vector(colMeans(class.x))
		
		sum_squares <- apply(class.x, 2, function(col) {
			(n_k - 1) * var(col)
		})
		
		list(xbar = xbar, sum_squares = sum_squares, n = n_k, pi_k = pi_k)
	})
	
	var_pool <- colSums(laply(estimators, function(class_est) class_est$sum_squares)) / N
	var_shrink <- var_shrinkage(N = N, K = num.classes, var_feature = var_pool, num_alphas = num_alphas, t = -1)
	
	estimators <- llply(estimators, function(class_estimators) {
		class_estimators$var <- var_shrink
		class_estimators
	})
	if(verbose) cat("Building SDLDA classifier...done!\n")
	
	sdlda.obj$N <- N
	sdlda.obj$classes <- levels(sdlda.obj$training$labels)
	sdlda.obj$estimators <- estimators
	
	class(sdlda.obj) <- "sdlda"
	
	sdlda.obj
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
