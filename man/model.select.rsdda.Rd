\name{model.select.rsdda}
\alias{model.select.rsdda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
TODO
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
model.select.rsdda(object, grid.size = 5, k = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{grid.size}{
%%     ~~Describe \code{grid.size} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
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
function(object, grid.size = 5, k = 1) {
	if (!inherits(object, "rsdda"))  {
		stop("object not of class 'rsdda'")
	}
	# Estimate lambda with Leave-k-out crossvalidation.
	if(k > 1) {
		k <- 1
		warning("Leave-k-out crossvalidation for k > 1 is not implemented yet. Using k = 1 instead...\n")
	}
	
	lambda.grid <- seq(0, 1, length = grid.size)
	
	# Here, we calculate the Leave-K-out Error Rates for each lambda in the grid.
	# NOTE: We are only simultaneously diagonalizing the covariance matrices ONCE.
	# If we did this for each observation left out using crossvalidation, the model
	# selection process would be take too long to run, especially in simulations
	# with a large number of replications.
	predictions <- laply(seq_len(object$N), function(i) {
		loo.rsdda.obj <- rsdda(object$training[-i, ], num.alphas = grid.size)
		laply(lambda.grid, function(lambda) {
			predict.rsdda(loo.rsdda.obj, object$training[i, -1], lambda = lambda)
		})
	})

	error.rates.loo <- apply(predictions, 2, function(predictions.lambda) {
		mean(object$training[, 1] != predictions.lambda)
	})
	
	# Find the lambdas that attain the minimum Leave-K-out error rate.
	optimal.lambdas <- lambda.grid[which(error.rates.loo == min(error.rates.loo))]
	
	# Finally, we add the optimal lambda to the RSDDA object.
	# If there is a tie for optimal lambda, randomly choose one from the optimal ones.
	object$lambda <- sample(optimal.lambdas, 1)
	
	return(object)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
