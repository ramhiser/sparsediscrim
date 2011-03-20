\name{risk_stein}
\alias{risk_stein}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
TODO
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
risk_stein(N, K, var_feature, num_alphas = 2, t = -1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
%%     ~~Describe \code{N} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{var_feature}{
%%     ~~Describe \code{var_feature} here~~
}
  \item{num_alphas}{
%%     ~~Describe \code{num_alphas} here~~
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
function(N, K, var_feature, num_alphas = 2, t = -1) {
	nu <- N - K
	p <- length(var_feature)
	alphas <- seq(0, 1, length = num_alphas)
	
	# The pooled variance is defined in Pang et al. (2009) as the geometric mean
	# of the sample variances of each feature.
	var_pool <- prod(var_feature)^(1 / p)
	
	# Here we compute the average risk for the Stein loss function on page 1023
	# for all values of alpha.
	risk_alphas <- laply(alphas, function(alpha) {
		risk <- h(nu = nu, p = p)^alpha * h(nu = nu, p = 1)^(1 - alpha)
		risk <- risk / (h(nu = nu, p = 1, t = alpha * t / p))^(p - 1)
		risk <- risk / h(nu = nu, p = 1, t = (1 - alpha + alpha / p) * t)
		risk <- risk * (var_pool)^(alpha * t)
		risk <- risk * mean(var_feature^(-alpha * t))
		risk <- risk - log(h(nu = nu, p = p)^alpha * h(nu = nu, p = 1)^(1 - alpha))
		risk <- risk - t * digamma(nu / 2)
		risk <- risk + t * log(nu / 2) - 1
	})

	# Which of the alphas empirically minimize this risk?
	# If there are ties in the minimum risk, we randomly select
	# the value of alpha from the minimizers.
	alpha.min.risk <- alphas[which(min(risk_alphas) == risk_alphas)]
	alpha.star <- sample(alpha.min.risk, 1)
	
	list(alpha = alpha.star, var_pool = var_pool)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
