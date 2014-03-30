# sparsediscrim

The R package `sparsediscrim` provides a collection of sparse and regularized discriminant
analysis classifiers that are especially useful for when applied to
small-sample, high-dimensional data sets.

The `sparsediscrim` package features the following classifier (the R function
is included within parentheses):

* High-Dimensional Regularized Discriminant Analysis (`hdrda`) from Ramey et al. (2014)

The `sparsediscrim` package also includes a variety of additional classifiers
intended for small-sample, high-dimensional data sets. These include:

* Diagonal Linear Discriminant Analysis (`dlda`) from Dudoit et al. (2002)
* Diagonal Quadratic Discriminant Analysis (`dqda`) from Dudoit et al. (2002)
* Linear Discriminant Analysis (LDA) with the Moore-Penrose Pseudo-Inverse (`lda_pseudo`)
* Linear Discriminant Analysis (LDA) with the Schafer-Strimmer estimator (`lda_schafer`)
* Linear Discriminant Analysis (LDA) with the Thomaz-Kitani-Gillies estimator (`lda_thomaz`)
* Minimum Distance Empirical Bayesian Estimator (`mdeb`) from Srivistava and Kubokawa (2007)
* Minimum Distance Rule using Modified Empirical Bayes (`mdmeb`) from Srivistava and Kubokawa (2007)
* Minimum Distance Rule using Moore-Penrose Inverse (`mdmp`) from Srivistava and Kubokawa (2007)
* Shrinkage-based Diagonal Linear Discriminant Analysis (`sdlda`) from Pang et al. (2009)
* Shrinkage-based Diagonal Quadratic Discriminant Analysis (`sdqda`) from Pang et al. (2009)
* Shrinkage-mean-based Diagonal Linear Discriminant Analysis (`smdlda`) from Tong et al. (2012)
* Shrinkage-mean-based Diagonal Quadratic Discriminant Analysis (`smdqda`) from Tong et al. (2012)

## Installation

You can install the stable version on [CRAN](http://cran.r-project.org/package=sparsediscrim):

```r
install.packages('sparsediscrim', dependencies = TRUE)
```

If you prefer to download the latest version, instead type:

```r
library(devtools)
install_github('sparsediscrim', 'ramhiser')
```
