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
sdlda(training.df, num.alphas = 5, jointdiag = "none", verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{training.df}{
%%     ~~Describe \code{training.df} here~~
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
function(training.df, num.alphas = 5, jointdiag = "none", verbose = FALSE, ...) {
	sdlda.obj <- list()
	sdlda.obj$training <- training.df
	
	if(jointdiag != "none") {
		if(verbose) cat("Simultaneously diagonalizing covariance matrices\n")
		joint.diag.out <- joint.diagonalization(sdlda.obj$training, method = jointdiag)
		sdlda.obj$training <- joint.diag.out$transformed.df
		sdlda.obj$jointdiag.B <- joint.diag.out$B
		sdlda.obj$jointdiag.method <- joint.diag.out$method
		if(verbose) cat("Simultaneously diagonalizing covariance matrices...done!\n")
	}
	
	if(verbose) cat("Building SDLDA classifier\n")
	N <- nrow(sdlda.obj$training)
	num.classes <- nlevels(sdlda.obj$training$labels)
	
	training.x <- as.matrix(sdlda.obj$training[,-1])
	dimnames(training.x) <- NULL
	
	estimators <- dlply(sdlda.obj$training, .(labels), function(class.df) {
		class.x <- as.matrix(class.df[,-1])
		dimnames(class.x) <- NULL
		
		n.k <- nrow(class.df)
		p.hat <- n.k / N
		xbar <- as.vector(colMeans(class.x))
		
		sum.squares <- apply(class.x, 2, function(col) {
			(n.k - 1) * var(col)
		})
		
		list(xbar = xbar, sum.squares = sum.squares, n = n.k, p.hat = p.hat)
	})
	
	var.pooled <- colSums(laply(estimators, function(class.est) class.est$sum.squares)) / N
	var.shrink <- var.shrinkage(N = N, K = num.classes, var.feature = var.pooled, num.alphas = num.alphas, t = -1)
	
	estimators <- llply(estimators, function(class.estimators) {
		class.estimators$var <- var.shrink
		class.estimators
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
