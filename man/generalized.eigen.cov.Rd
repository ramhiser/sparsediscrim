\name{generalized.eigen.cov}
\alias{generalized.eigen.cov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generalized.Eigen.Cov
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
generalized.eigen.cov(df, shrink = FALSE, shrink.val = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
%%     ~~Describe \code{df} here~~
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
function(df, shrink = FALSE, shrink.val = 0.01) {
K <- nlevels(factor(df$labels))
stopifnot(K == 2)
if(shrink == FALSE) {
shrink.val <- 0
  }
covs <- dlply(df, .(labels), function(df_k) {
n_k <- nrow(df_k)
p <- ncol(df_k) - 1
cov.k <- (1 - shrink.val) * (n_k - 1) * cov(df_k[,-1]) / n_k + shrink.val * diag(p)
dimnames(cov.k) <- NULL
cov.k
  })
covs
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
