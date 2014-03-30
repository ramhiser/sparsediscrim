# sparsediscrim

The R package `sparsediscrim` provides a collection of sparse and regularized discriminant
analysis classifiers that are especially useful for when applied to
small-sample, high-dimensional data sets.

The `sparsediscrim` package features the following classifier (the R function
is included within parentheses):

* High-Dimensional Regularized Discriminant Analysis (`hdrda`) from Ramey et al. (2014)

The `sparsediscrim` package also includes a variety of additional classifiers
intended for small-sample, high-dimensional data sets. These include:

* Diagonal Linear Discriminant Analysis from Dudoit et al. (2002) (`dlda`)
* Diagonal Quadratic Discriminant Analysis from Dudoit et al. (2002) (`dqda`)
* Linear Discriminant Analysis (LDA) with the Moore-Penrose Pseudo-Inverse (`lda_pseudo`)
* Linear Discriminant Analysis (LDA) with the Schafer-Strimmer estimator (`lda_schafer`)
* Linear Discriminant Analysis (LDA) with the Thomaz-Kitani-Gillies estimator (`lda_thomaz`)
* Minimum Distance Empirical Bayesian Estimator from Srivistava and Kubokawa (2007) (`mdeb`)
* Minimum Distance Rule using Modified Empirical Bayes from Srivistava and Kubokawa (2007) (`mdmeb`)
* Minimum Distance Rule using Moore-Penrose Inverse from Srivistava and Kubokawa (2007) (`mdmp`)
* Shrinkage-based Diagonal Linear Discriminant Analysis from Pang et al. (2009) (`sdlda`)
* Shrinkage-based Diagonal Quadratic Discriminant Analysis from Pang et al. (2009) (`sdqda`)
* Shrinkage-mean-based Diagonal Linear Discriminant Analysis from Tong et al. (2012) (`smdlda`)
* Shrinkage-mean-based Diagonal Quadratic Discriminant Analysis from Tong et al. (2012) (`smdqda`)

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
