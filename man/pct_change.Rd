\name{pct_change}
\alias{pct_change}
\title{Calculates the percentage change for each element in the given vector}
\usage{
  pct_change(x)
}
\arguments{
  \item{x}{vector}
}
\value{
  vector of length \code{length(x) - 1}
}
\description{
  For each element in the vector, \code{x}, we calculate
  the percentage change.
}
\details{
  The returned vector has one less observation than
  \code{x}.
}
\examples{
x <- seq_len(10)^2
out <- pct_change(x)
out
length(out) == 9
}

