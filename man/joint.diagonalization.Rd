\name{joint.diagonalization}
\alias{joint.diagonalization}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
TODO
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
joint.diagonalization(df, method = c("none", "general-eigen", "asfari", "jade", "jedi", "qdiag", "ffdiag", "jadiag", "uwedge"), tol = 1e-06, max.iter = 250, shrink = TRUE, shrink.val = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
%%     ~~Describe \code{df} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{tol}{
%%     ~~Describe \code{tol} here~~
}
  \item{max.iter}{
%%     ~~Describe \code{max.iter} here~~
}
  \item{shrink}{
%%     ~~Describe \code{shrink} here~~
}
  \item{shrink.val}{
%%     ~~Describe \code{shrink.val} here~~
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
function(df, method = c("none", "general-eigen", "asfari", "jade", "jedi", "qdiag", "ffdiag", "jadiag", "uwedge"),
	tol = 1e-6, max.iter = 250, shrink = TRUE, shrink.val = 0.01) {

	method <- match.arg(method)
	if(method == "none") {
		B <- diag(ncol(df) - 1)
	} else if(method == "general-eigen") {
		warning("Generalized Eigenvalues has not been implemented yet.")
	} else if(method == "asfari") {
		asfari.cov <- asfari.cov(df, shrink = TRUE, shrink.val = shrink.val)
		B <- LUJID(asfari.cov, mode = 'B', ERR = tol, RBALANCE = 3, ITER = max.iter)
		df <- joint.diagonalization.transform(df, B)
	} else if(method == "jade") {
		warning("JADE algorithm has not been implemented yet.")
		B <- diag(ncol(df) - 1)
	} else if(method == "jedi") {
		ajd.cov <- jointDiag.cov(df, shrink = TRUE, shrink.val = shrink.val)
		B <- jedi(ajd.cov, eps = tol, itermax = max.iter)$A
		df <- joint.diagonalization.transform(df, B)
	} else {
		ajd.cov <- jointDiag.cov(df, shrink = TRUE, shrink.val = shrink.val)
		B <- ajd(ajd.cov, method = method, eps = tol, itermax = max.iter)$B
		df <- joint.diagonalization.transform(df, B)
	}
	list(transformed.df = df, B = B, method = method)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
