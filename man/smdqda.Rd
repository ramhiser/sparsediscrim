\name{smdqda}
\alias{predict.smdqda}
\alias{print.smdqda}
\alias{smdqda}
\alias{smdqda.default}
\alias{smdqda.formula}
\title{Shrinkage-mean-based Diagonal Quadratic Discriminant Analysis (SmDQDA) from
Tong, Chen, and Zhao (2012)}
\usage{
smdqda(x, ...)

\method{smdqda}{default}(x, y, prior = NULL, ...)

\method{smdqda}{formula}(formula, data, prior = NULL, ...)

\method{print}{smdqda}(x, ...)

\method{predict}{smdqda}(object, newdata, ...)
}
\arguments{
  \item{x}{matrix containing the training data. The rows
  are the sample observations, and the columns are the
  features.}

  \item{y}{vector of class labels for each training
  observation}

  \item{prior}{vector with prior probabilities for each
  class. If NULL (default), then equal probabilities are
  used. See details.}

  \item{formula}{A formula of the form \code{groups ~ x1 +
  x2 + ...} That is, the response is the grouping factor
  and the right hand side specifies the (non-factor)
  discriminators.}

  \item{data}{data frame from which variables specified in
  \code{formula} are preferentially to be taken.}

  \item{x}{object to print}

  \item{...}{unused}

  \item{object}{trained SmDQDA object}

  \item{newdata}{matrix of observations to predict. Each
  row corresponds to a new observation.}

  \item{...}{additional arguments}
}
\value{
\code{smdqda} object that contains the trained SmDQDA
classifier

list predicted class memberships of each row in newdata
}
\description{
Given a set of training data, this function builds the
Shrinkage-mean-based Diagonal Quadratic Discriminant
Analysis (SmDQDA) classifier from Tong, Chen, and Zhao
(2012). The SmDQDA classifier incorporates a Lindley-type
shrunken mean estimator into the DQDA classifier from
Dudoit et al. (2002). For more about the DQDA classifier,
see \code{\link{dqda}}.

Summarizes the trained SmDQDA classifier in a nice manner.

The SmDQDA classifier is a modification to QDA, where the
off-diagonal elements of the pooled sample covariance
matrix are set to zero.
}
\details{
The DQDA classifier is a modification to the well-known QDA
classifier, where the off-diagonal elements of each class
covariance matrix are assumed to be zero -- the features
are assumed to be uncorrelated. Under multivariate
normality, the assumption uncorrelated features is
equivalent to the assumption of independent features. The
feature-independence assumption is a notable attribute of
the Naive Bayes classifier family. The benefit of these
classifiers is that they are fast and have much fewer
parameters to estimate, especially when the number of
features is quite large.

The matrix of training observations are given in \code{x}.
The rows of \code{x} contain the sample observations, and
the columns contain the features for each training
observation.

The vector of class labels given in \code{y} are coerced to
a \code{factor}. The length of \code{y} should match the
number of rows in \code{x}.

An error is thrown if a given class has less than 2
observations because the variance for each feature within a
class cannot be estimated with less than 2 observations.

The vector, \code{prior}, contains the \emph{a priori}
class membership for each class. If \code{prior} is NULL
(default), the class membership probabilities are estimated
as the sample proportion of observations belonging to each
class. Otherwise, \code{prior} should be a vector with the
same length as the number of classes in \code{y}. The
\code{prior} probabilties should be nonnegative and sum to
one.
}
\examples{
n <- nrow(iris)
train <- sample(seq_len(n), n / 2)
smdqda_out <- smdqda(Species ~ ., data = iris[train, ])
predicted <- predict(smdqda_out, iris[-train, -5])$class

smdqda_out2 <- smdqda(x = iris[train, -5], y = iris[train, 5])
predicted2 <- predict(smdqda_out2, iris[-train, -5])$class
all.equal(predicted, predicted2)
}
\references{
Tong, T., Chen, L., and Zhao, H. (2012), "Improved Mean
Estimation and Its Application to Diagonal Discriminant
Analysis," Bioinformatics, 28, 4, 531-537.
\url{http://bioinformatics.oxfordjournals.org/content/28/4/531.long}

Dudoit, S., Fridlyand, J., & Speed, T. P. (2002).
"Comparison of Discrimination Methods for the
Classification of Tumors Using Gene Expression Data,"
Journal of the American Statistical Association, 97, 457,
77-87.

Dudoit, S., Fridlyand, J., & Speed, T. P. (2002).
"Comparison of Discrimination Methods for the
Classification of Tumors Using Gene Expression Data,"
Journal of the American Statistical Association, 97, 457,
77-87.
}

