\name{predict.dqda}
\alias{predict.dqda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
TODO
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
predict.dqda(object, newdata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{newdata}{
%%     ~~Describe \code{newdata} here~~
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
function(object, newdata) {
	if (!inherits(object, "dqda"))  {
		stop("object not of class 'dqda'")
	}
	newdata <- as.matrix(newdata)
	dimnames(newdata) <- NULL
	
	predictions <- apply(newdata, 1, function(obs) {
		if(!is.null(object$jointdiag.method) && object$jointdiag.method != "none") {
			obs <- obs \%*\% t(object$jointdiag.B)
		}
		scores <- sapply(object$estimators, function(class.est) {
			sum((obs - class.est$xbar)^2 / class.est$var) + sum(log(class.est$var)) - 2 * log(class.est$p.hat)
		})
		predicted.class <- object$classes[which.min(scores)]
		predicted.class
	})
	
	predictions
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
