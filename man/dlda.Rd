\name{dlda}
\alias{dlda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
TODO
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dlda(training.df, jointdiag = "none", verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{training.df}{
%%     ~~Describe \code{training.df} here~~
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
function(training.df, jointdiag = "none", verbose = FALSE, ...) {
	dlda.obj <- list()
	dlda.obj$training <- training.df
	
	if(jointdiag != "none") {
		if(verbose) message("Simultaneously diagonalizing covariance matrices... ", appendLF = FALSE)
		joint.diag.out <- joint.diagonalization(dlda.obj$training, method = jointdiag)
		dlda.obj$training <- joint.diag.out$transformed.df
		dlda.obj$jointdiag.B <- joint.diag.out$B
		dlda.obj$jointdiag.method <- joint.diag.out$method
		if(verbose) message("done!")
	}

	if(verbose) cat("Building DLDA classifier\n")
	N <- nrow(dlda.obj$training)
	
	training.x <- as.matrix(dlda.obj$training[,-1])
	dimnames(training.x) <- NULL
	
	estimators <- dlply(dlda.obj$training, .(labels), function(class.df) {
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
	
	estimators <- llply(estimators, function(class.estimators) {
		class.estimators$var <- var.pooled
		class.estimators
	})
	if(verbose) cat("Building DLDA classifier...done!\n")
	
	dlda.obj$N <- N
	dlda.obj$classes <- levels(dlda.obj$training$labels)
	dlda.obj$estimators <- estimators
	
	class(dlda.obj) <- "dlda"
	
	dlda.obj
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
