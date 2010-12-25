\name{var.shrinkage}
\alias{var.shrinkage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
TODO
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
var.shrinkage(N, K, var.feature, num.alphas = 2, t = -1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
%%     ~~Describe \code{N} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{var.feature}{
%%     ~~Describe \code{var.feature} here~~
}
  \item{num.alphas}{
%%     ~~Describe \code{num.alphas} here~~
}
  \item{t}{
%%     ~~Describe \code{t} here~~
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
function(N, K, var.feature, num.alphas = 2, t = -1) {
	nu <- N - K
	p <- length(var.feature)
	
	risk.stein.out <- risk.stein(N = N, K = K, var.feature = var.feature, num.alphas = num.alphas, t = t)

	var.pooled <- risk.stein.out$var.pooled
	alpha <- risk.stein.out$alpha
	
	var.feature.shrink <- (h(nu = nu, p = p, t = t) * var.pooled)^alpha * (h(nu = nu, p = 1, t = t) * var.feature)^(1 - alpha)
	var.feature.shrink
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
