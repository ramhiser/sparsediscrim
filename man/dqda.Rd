\name{dqda}
\alias{dqda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
TODO
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dqda(train_df, jointdiag = "none", verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{train_df}{
%%     ~~Describe \code{train_df} here~~
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
function(train_df, jointdiag = "none", verbose = FALSE, ...) {
	dqda.obj <- list()
	dqda.obj$training <- train_df
	
	if(jointdiag != "none") {
		if(verbose) message("Simultaneously diagonalizing covariance matrices... ", appendLF = FALSE)
		joint.diag.out <- joint.diagonalization(dqda.obj$training, method = jointdiag)
		dqda.obj$training <- joint.diag.out$transformed.df
		dqda.obj$jointdiag.B <- joint.diag.out$B
		dqda.obj$jointdiag.method <- joint.diag.out$method
		if(verbose) message("done!")
	}
	
	if(verbose) cat("Building DQDA classifier\n")
	N <- nrow(dqda.obj$training)
	
	estimators <- dlply(dqda.obj$training, .(labels), function(df_k) {
		n_k <- nrow(df_k)
		pi_k <- n_k / N
		xbar <- as.vector(colMeans(df_k[, -1]))
		var <- apply(df_k[,-1], 2, function(col) {
			(n_k - 1) * var(col) / n_k
		})
		list(xbar = xbar, var = var, n = n_k, pi_k = pi_k)
	})
	
	if(verbose) cat("Building DQDA classifier...done!\n")
	
	dqda.obj$N <- N
	dqda.obj$classes <- levels(dqda.obj$training$labels)
	dqda.obj$estimators <- estimators
	
	class(dqda.obj) <- "dqda"
	
	dqda.obj
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
