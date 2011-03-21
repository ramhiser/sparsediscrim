\name{predict.rsdda}
\alias{predict.rsdda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
TODO
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
predict.rsdda(object, newdata, num_lambdas = 5, lambda = NULL, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{newdata}{
%%     ~~Describe \code{newdata} here~~
}
  \item{num_lambdas}{
%%     ~~Describe \code{num_lambdas} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
function(object, newdata, num_lambdas = 5, lambda = NULL, verbose = FALSE) {
	if (!inherits(object, "rsdda"))  {
		stop("object not of class 'rsdda'")
	}
	
	if(is.null(lambda) || !is.numeric(lambda)) {
		if(verbose) cat("Model selection for RSDDA\n")
		object <- model.select.rsdda(object, grid_size = num_lambdas)
		if(verbose) cat("Model selection for RSDDA...done!\n")
		if(verbose) cat("RSDDA Lambda Parameter:", object$lambda, "\n")
	}
	else {
		object$lambda <- lambda
	}
	
	newdata <- as.matrix(newdata)
	dimnames(newdata) <- NULL
	
	predictions <- apply(newdata, 1, function(obs) {
		if(!is.null(object$jointdiag.method) && object$jointdiag.method != "none") {
			obs <- obs \%*\% t(object$jointdiag.B)
		}
		scores <- sapply(object$estimators, function(class_est) {
			var.rsdda <- (class_est$var.k)^(1 - object$lambda) * (class_est$var.pool)^(object$lambda)
			sum((obs - class_est$xbar)^2 * var.rsdda) - sum(log(var.rsdda)) - 2 * log(class_est$pi_k)
		})
		
		prediction <- object$classes[which.min(scores)]
		prediction
	})

	predictions
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
