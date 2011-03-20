\name{sdqda}
\alias{sdqda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
TODO
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sdqda(train_df, num_alphas = 5, jointdiag = "none", verbose = FALSE, ...)
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
	sdqda.obj <- list()
	sdqda.obj$training <- train_df
	
	if(jointdiag != "none") {
		if(verbose) message("Simultaneously diagonalizing covariance matrices... ", appendLF = FALSE)
		joint.diag.out <- joint.diagonalization(sdqda.obj$training, method = jointdiag)
		sdqda.obj$training <- joint.diag.out$transformed.df
		sdqda.obj$jointdiag.B <- joint.diag.out$B
		sdqda.obj$jointdiag.method <- joint.diag.out$method
		if(verbose) message("done!")
	}
	
	if(verbose) cat("Building SDQDA classifier\n")
	N <- nrow(sdqda.obj$training)
	num.classes <- nlevels(sdqda.obj$training$labels)
	
	estimators <- dlply(sdqda.obj$training, .(labels), function(df_k) {
		class.x <- as.matrix(df_k[, -1])
		dimnames(class.x) <- NULL
		
		n_k <- nrow(class.x)
		pi_k <- n_k / N
		
		xbar <- as.vector(colMeans(class.x))
		var <- apply(class.x, 2, function(col) {
			(n_k - 1) * var(col) / n_k
		})
		
		var_shrink <- var_shrinkage(N = n_k, K = 1, var_feature = var, num_alphas = num_alphas, t = -1)
		
		list(xbar = xbar, var = var_shrink, n = n_k, pi_k = pi_k)
	})
	if(verbose) cat("Building SDQDA classifier...done!\n")
	
	sdqda.obj$N <- N
	sdqda.obj$classes <- levels(sdqda.obj$training$labels)
	sdqda.obj$estimators <- estimators
	
	class(sdqda.obj) <- "sdqda"
	
	sdqda.obj
	
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
