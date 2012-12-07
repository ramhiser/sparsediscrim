\name{yj_inverse}
\alias{yj_inverse}
\title{Computes the inverse of Yeo-Johnson transformed random variate(s)}
\usage{
  yj_inverse(y, lambda, tol = sqrt(.Machine$double.eps))
}
\arguments{
  \item{y}{vector of random variate(s)}

  \item{lambda}{numeric Yeo-Johnson power transformation
  parameter value}

  \item{tol}{numeric tolerance value that determines
  numerical equality}
}
\value{
  the inverse of the Yeo-Johnson random variate, \code{y}
}
\description{
  TODO: Briefly describe Yeo-Johnson is generalization of
  Box-Cox TODO: Provide the formulas for Yeo-Johnson
  transformation (given in paper) TODO: Provide the inverse
  formulas for Yeo-Johnson (on whiteboard and iPhone) TODO:
  Cite Yeo-Johnson paper TODO: Reference the 'car' package
  to find the fitted 'lambda'
}

